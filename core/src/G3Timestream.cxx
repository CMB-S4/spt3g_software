#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <G3Timestream.h>
#include <G3Units.h>

#ifdef G3_HAS_FLAC
#include <FLAC/stream_encoder.h>
#include <cmath>

template<typename A>
struct FlacDecoderCallbackArgs {
	A *inbuf;
	std::vector<int32_t> *outbuf;
	size_t pos;
	size_t nbytes;
};

enum FLACNaNFlag {
	NoNan = 0,
	AllNan = 1,
	SomeNan = 2
};

static FLAC__StreamEncoderWriteStatus flac_encoder_write_cb(
    const FLAC__StreamEncoder *encoder, const FLAC__byte buffer[], size_t bytes,
    unsigned samples, unsigned current_frame, void *client_data)
{
	std::vector<uint8_t> *outbuf = (std::vector<uint8_t> *)(client_data);

	outbuf->insert(outbuf->end(), buffer, buffer + bytes);
	return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

template<typename A>
static FLAC__StreamDecoderReadStatus flac_decoder_read_cb(
    const FLAC__StreamDecoder *decoder, FLAC__byte buffer[], size_t *bytes,
    void *client_data)
{
	FlacDecoderCallbackArgs<A> *args =
	    (FlacDecoderCallbackArgs<A> *)(client_data);

	ssize_t bytes_left = ssize_t(args->nbytes) - args->pos;

	if (bytes_left <= 0 || *bytes == 0) {
		*bytes = 0;
		return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM;
	} else if (*bytes >= size_t(bytes_left)) {
		*bytes = bytes_left;
		args->inbuf->template loadBinary<1>(buffer, bytes_left);
		args->pos += bytes_left;
		return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
	} else {
		args->inbuf->template loadBinary<1>(buffer, *bytes);
		args->pos += *bytes;
		return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
	}
}
  
template<typename A>
static FLAC__StreamDecoderWriteStatus flac_decoder_write_cb(
    const FLAC__StreamDecoder *decoder, const FLAC__Frame *frame,
    const FLAC__int32 *const buffer[], void *client_data)
{
	FlacDecoderCallbackArgs<A> *args =
	    (FlacDecoderCallbackArgs<A> *)(client_data);

	size_t oldsize = args->outbuf->size();
	args->outbuf->resize(oldsize + frame->header.blocksize);
	for (size_t i = 0; i < frame->header.blocksize; i++)
		(*args->outbuf)[oldsize + i] = buffer[0][i];
	return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}
  
static void flac_decoder_error_cb(const FLAC__StreamDecoder *decoder,
    FLAC__StreamDecoderErrorStatus status, void *client_data)
{

	switch (status) {
	case FLAC__STREAM_DECODER_ERROR_STATUS_LOST_SYNC:
		log_fatal("FLAC decoding error (lost sync)");
	case FLAC__STREAM_DECODER_ERROR_STATUS_BAD_HEADER:
		log_fatal("FLAC decoding error (bad header)");
	case FLAC__STREAM_DECODER_ERROR_STATUS_FRAME_CRC_MISMATCH:
		log_fatal("FLAC decoding error (CRC mismatch)");
	case FLAC__STREAM_DECODER_ERROR_STATUS_UNPARSEABLE_STREAM:
		log_fatal("FLAC decoding error (unparseable stream)");
	default:
		log_fatal("FLAC decoding error (%d)", status);
	}
}

extern "C"{
	// Provide our own declaration of this function.
	// This libFLAC interface is private but stable, and this use is officially sanctioned:
	// https://github.com/xiph/flac/commit/3baaf23faa05eca1cfc34737d95131ad0b628d4c
	FLAC__bool FLAC__stream_encoder_set_do_md5(FLAC__StreamEncoder *encoder, FLAC__bool value);
}
#endif

template <class A> void G3Timestream::save(A &ar, unsigned v) const
{
	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("units", units);
	ar & cereal::make_nvp("start", start);
	ar & cereal::make_nvp("stop", stop);
	ar & cereal::make_nvp("flac", use_flac_);

#ifdef G3_HAS_FLAC
	if (use_flac_) {
		std::vector<int32_t> inbuf;
		std::vector<uint8_t> outbuf;
		const int32_t *chanmap[1];
		uint8_t nanflag;
		size_t nans = 0;

		if (units != Counts && units != None)
			log_fatal("Cannot use FLAC on non-counts timestreams");

		DataType data_type_out = data_type_;

		// Copy to 24-bit integers
		inbuf.resize(size());
		switch (flac_depth_) {
		case 24:
			switch (data_type_) {
			case TS_DOUBLE:
				data_type_out = TS_FLOAT;
				for (size_t i = 0; i < size(); i++)
					inbuf[i] = ((int32_t(((double *)data_)[i]) & 0x00ffffff) << 8) >> 8;
				break;
			case TS_FLOAT:
				for (size_t i = 0; i < size(); i++)
					inbuf[i] = ((int32_t(((float *)data_)[i]) & 0x00ffffff) << 8) >> 8;
				break;
			case TS_INT32:
				{
					// Using this rather raw form for the loop can enable automatic
					// unrolling and vectorization.
					int32_t* in_ptr=(int32_t *)data_;
					int32_t* out_ptr=&inbuf[0];
					for(int32_t* end=in_ptr+size(); in_ptr!=end; in_ptr++,out_ptr++)
						*out_ptr = ((*in_ptr & 0x00ffffff) << 8) >> 8;
				}
				break;
			case TS_INT64:
				data_type_out = TS_INT32;
				for (size_t i = 0; i < size(); i++)
					inbuf[i] = ((int32_t(((int64_t *)data_)[i]) & 0x00ffffff) << 8) >> 8;
				break;
			default:
				log_fatal("Unknown timestream datatype %d", data_type_);
			}
			break;
		case 32:
			if (FLAC_API_VERSION_CURRENT < 13)
				log_fatal("32-bit compression is not supported, "
				    "please upgrade FLAC to version 1.4 or newer");

			switch (data_type_) {
			case TS_DOUBLE:
				for (size_t i = 0; i < size(); i++)
					inbuf[i] = int32_t(((double *)data_)[i]);
				break;
			case TS_FLOAT:
				for (size_t i = 0; i < size(); i++)
					inbuf[i] = int32_t(((float *)data_)[i]);
				break;
			case TS_INT32:
				memcpy(inbuf.data(), data_, size() * sizeof(int32_t));
				break;
			case TS_INT64:
				for (size_t i = 0; i < size(); i++)
					inbuf[i] = int32_t(((int64_t *)data_)[i]);
				break;
			default:
				log_fatal("Unknown timestream datatype %d", data_type_);
			}
			break;
		default:
			log_fatal("Invalid FLAC bit depth %d", flac_depth_);
		}
		chanmap[0] = &inbuf[0];

		ar & cereal::make_nvp("flac_depth", flac_depth_);
		ar & cereal::make_nvp("data_type", data_type_out);

		// Mark bad samples using an out-of-band signal. Since we
		// convert to 24-bit integers going into FLAC, which have no
		// out-of-range values for signalling, this requires special
		// care. Usually a timestream is either all-valid or
		// all-invalid, so signal that with a single byte flag. In the
		// rare case that only some samples are valid, store a
		// validity mask.
		std::vector<bool> nanbuf(size(), false);
		if(data_type_==TS_DOUBLE || data_type_==TS_FLOAT){
			for (size_t i = 0; i < size(); i++) {
				if (!std::isfinite((*this)[i])) {
					nans++;
					nanbuf[i] = true;
					inbuf[i] = 0;
				}
			}
		}
		nanflag = SomeNan;
		if (nans == 0)
			nanflag = NoNan;
		else if (nans == size())
			nanflag = AllNan;
		ar & cereal::make_nvp("nanflag", nanflag);
		if (nanflag == SomeNan)
			ar & cereal::make_nvp("nanmask", nanbuf);

		// Now do FLAC encoding
		FLAC__StreamEncoder *encoder = FLAC__stream_encoder_new();
		FLAC__stream_encoder_set_channels(encoder, 1);
		FLAC__stream_encoder_set_bits_per_sample(encoder, flac_depth_);
		FLAC__stream_encoder_set_compression_level(encoder, use_flac_);
		FLAC__stream_encoder_set_do_md5(encoder, false);
		FLAC__stream_encoder_init_stream(encoder,
		    flac_encoder_write_cb, NULL, NULL, NULL, (void*)(&outbuf));
		FLAC__stream_encoder_process (encoder, chanmap, inbuf.size());
		FLAC__stream_encoder_finish(encoder);
		FLAC__stream_encoder_delete(encoder);

		ar & cereal::make_nvp("data", outbuf);
	} else {
#endif
		ar & cereal::make_nvp("data_type", data_type_);
		if (buffer_) {
			ar & cereal::make_nvp("data", *buffer_);
		} else {
			switch (data_type_) {
			case TS_DOUBLE: {
				std::vector<double> data((double *)data_,
				    (double *)data_ + len_);
				ar & cereal::make_nvp("data", data);
				break;
			}
			case TS_FLOAT: {
				std::vector<float> data((float *)data_,
				    (float *)data_ + len_);
				ar & cereal::make_nvp("data", data);
				break;
			}
			case TS_INT32: {
				std::vector<int32_t> data((int32_t *)data_,
				    (int32_t *)data_ + len_);
				ar & cereal::make_nvp("data", data);
				break;
			}
			case TS_INT64: {
				std::vector<int64_t> data((int64_t *)data_,
				    (int64_t *)data_ + len_);
				ar & cereal::make_nvp("data", data);
				break;
			}
			default:
				log_fatal("Unknown timestream datatype %d", data_type_);
			}
		}
#ifdef G3_HAS_FLAC
	}
#endif
}

template <typename T>
std::vector<T> *
unpack_flac(const std::vector<int32_t> &buf, uint8_t nanflag, const std::vector<bool> &nanbuf)
{
	// Represent data as floats internally if possible, to allow NaNs,
	// which we try to pull through the process to signal missing data.

	// Convert data format
	auto data = new std::vector<T>(buf.size());
	for (size_t i = 0; i < buf.size(); i++)
		(*data)[i] = buf[i];

	// Apply NaN mask
	if (nanflag == AllNan) {
		for (size_t i = 0; i < buf.size(); i++)
			(*data)[i] = NAN;
	} else if (nanflag == SomeNan) {
		for (size_t i = 0; i < buf.size(); i++)
			if (nanbuf[i])
				(*data)[i] = NAN;
	}

	return data;
}

template <class A> void G3Timestream::load(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("units", units);
	if (v >= 2) {
		ar & cereal::make_nvp("start", start);
		ar & cereal::make_nvp("stop", stop);
	}
	ar & cereal::make_nvp("flac", use_flac_);

	if (use_flac_) {
#ifdef G3_HAS_FLAC
		FlacDecoderCallbackArgs<A> callback;
		uint8_t nanflag;
		std::vector<bool> nanbuf;

		callback.inbuf = &ar;
		if (buffer_)
			delete buffer_;
		root_data_ref_.reset();
		buffer_ = NULL;
		callback.outbuf = new std::vector<int32_t>();
		callback.pos = 0;

		if (units != Counts && units != None)
			log_fatal("Cannot use FLAC on non-counts timestreams");

		if (v >= 4) {
			ar & cereal::make_nvp("flac_depth", flac_depth_);
			if (flac_depth_ > 24 && FLAC_API_VERSION_CURRENT < 13)
				log_fatal("32-bit decompression is not supported, "
				    "please upgrade FLAC to version 1.4 or newer");
			ar & cereal::make_nvp("data_type", data_type_);
		} else {
			flac_depth_ = 24;
			data_type_ = TS_FLOAT;
		}

		ar & cereal::make_nvp("nanflag", nanflag);
		if (nanflag == SomeNan)
			ar & cereal::make_nvp("nanmask", nanbuf);

		ar & cereal::make_size_tag(callback.nbytes);

		// Typical compression ratio: N bytes in input = N samples
		callback.outbuf->reserve(callback.nbytes);

		FLAC__StreamDecoder *decoder = FLAC__stream_decoder_new();
		FLAC__stream_decoder_set_md5_checking(decoder, false);
		FLAC__stream_decoder_init_stream(decoder,
		    flac_decoder_read_cb<A>, NULL, NULL, NULL, NULL,
		    flac_decoder_write_cb<A>, NULL, flac_decoder_error_cb,
		    (void*)(&callback));
		FLAC__stream_decoder_process_until_end_of_stream(decoder);
		FLAC__stream_decoder_finish(decoder);
		FLAC__stream_decoder_delete(decoder);

		len_ = callback.outbuf->size();

		switch (data_type_) {
		case TS_INT32:
			// Short-circuit for int32
			root_data_ref_ = std::shared_ptr<std::vector<int32_t> >(
			    callback.outbuf);
			data_ = &(*callback.outbuf)[0];
			return;
		case TS_INT64: {
			auto data = new std::vector<int64_t>(size());
			for (size_t i = 0; i < size(); i++)
				(*data)[i] = (*callback.outbuf)[i];
			root_data_ref_ = std::shared_ptr<std::vector<int64_t> >(data);
			data_ = &(*data)[0];
			break;
		}
		case TS_FLOAT: {
			auto data = unpack_flac<float>(*callback.outbuf, nanflag, nanbuf);
			root_data_ref_ = std::shared_ptr<std::vector<float> >(data);
			data_ = &(*data)[0];
			break;
		}
		case TS_DOUBLE:
			buffer_ = unpack_flac<double>(*callback.outbuf, nanflag, nanbuf);
			data_ = &(*buffer_)[0];
			break;
		default:
			log_fatal("Unknown timestream datatype %d", data_type_);
		}

		delete callback.outbuf;
#else
		log_fatal("Trying to read FLAC-compressed timestreams but built without FLAC support");
#endif
	} else {
		if (buffer_)
			delete buffer_;
		buffer_ = NULL;
		root_data_ref_.reset();

		if (v >= 3)
			ar & cereal::make_nvp("data_type", data_type_);
		else
			data_type_ = TS_DOUBLE;
		switch (data_type_) {
		case TS_DOUBLE:
			buffer_ = new std::vector<double>();
			ar & cereal::make_nvp("data", *buffer_);
			len_ = buffer_->size();
			data_ = &(*buffer_)[0];
			break;
		case TS_FLOAT: {
			std::vector<float> *data = new std::vector<float>();
			ar & cereal::make_nvp("data", *data);
			root_data_ref_ = std::shared_ptr<std::vector<float> >(
			    data);
			len_ = data->size();
			data_ = &(*data)[0];
			break;
			}
		case TS_INT32: {
			std::vector<int32_t> *data = new std::vector<int32_t>();
			ar & cereal::make_nvp("data", *data);
			root_data_ref_ = std::shared_ptr<
			    std::vector<int32_t> >(data);
			len_ = data->size();
			data_ = &(*data)[0];
			break;
		}
		case TS_INT64: {
			std::vector<int64_t> *data = new std::vector<int64_t>();
			ar & cereal::make_nvp("data", *data);
			root_data_ref_ = std::shared_ptr<
			    std::vector<int64_t> >(data);
			len_ = data->size();
			data_ = &(*data)[0];
			break;
		}
		default:
			log_fatal("Unknown timestream datatype %d", data_type_);
		}
	}
}

G3Timestream::G3Timestream(const G3Timestream &r) :
    units(r.units), start(r.start), stop(r.stop), use_flac_(r.use_flac_),
    flac_depth_(r.flac_depth_), len_(r.len_), data_type_(r.data_type_)
{
	// Copy constructor needs to copy data, which always involves
	// allocating the internal buffer.
	if (r.buffer_) {
		buffer_ = new std::vector<double>(*r.buffer_);
		data_ = (&(*buffer_)[0]);
	} else if (r.data_type_ == TS_DOUBLE) {
		buffer_ = new std::vector<double>(len_);
		for (size_t i = 0; i < len_; i++)
			(*buffer_)[i] = r[i];
		data_ = (&(*buffer_)[0]);
	} else {
		buffer_ = NULL;
		size_t element = 0;
		switch (data_type_) {
		case TS_DOUBLE:
			__builtin_unreachable(); // Handled above
		case TS_FLOAT: {
			std::vector<float> *data = new std::vector<float>(len_);
			root_data_ref_ = std::shared_ptr<std::vector<float> >(
			    data);
			data_ = &(*data)[0];
			element = 4;
			break;
			}
		case TS_INT32: {
			std::vector<int32_t> *data =
			    new std::vector<int32_t>(len_);
			root_data_ref_ = std::shared_ptr<
			    std::vector<int32_t> >(data);
			data_ = &(*data)[0];
			element = 4;
			break;
		}
		case TS_INT64: {
			std::vector<int64_t> *data =
			    new std::vector<int64_t>(len_);
			root_data_ref_ = std::shared_ptr<
			    std::vector<int64_t> >(data);
			data_ = &(*data)[0];
			element = 8;
			break;
		}
		default:
			log_fatal("Unknown timestream datatype %d", data_type_);
		}
		memcpy(data_, r.data_, element*len_);
	}
}

double G3Timestream::GetSampleRate() const
{
	return double(size() - 1)/double(stop.time - start.time);
}

void G3Timestream::SetFLACCompression(int use_flac)
{

#ifdef G3_HAS_FLAC
	if (use_flac != 0 && units != Counts && units != None)
		log_fatal("Cannot use FLAC on non-counts timestreams");

	use_flac_ = use_flac;
#else
	if (use_flac != 0)
		log_fatal("Built without FLAC support");
#endif
}

void G3Timestream::SetFLACBitDepth(int bit_depth)
{
	if (bit_depth != 24 && bit_depth != 32)
		log_fatal("Invalid flac bit depth %d", bit_depth);

	flac_depth_ = bit_depth;
}

G3Timestream G3Timestream::operator +(const G3Timestream &r) const
{
	G3Timestream ret(*this);

	if (r.size() != size())
		log_fatal("Adding timestreams of unequal length");
	if (r.units != units && r.units != None && units != None)
		log_fatal("Adding timestreams of unequal units");
	for (size_t i = 0; i < size(); i++)
		ret[i] = (*this)[i] + r[i];

	return ret;
}

G3Timestream &G3Timestream::operator +=(const G3Timestream &r)
{
	if (r.size() != size())
		log_fatal("Adding timestreams of unequal length");
	if (r.units != units && r.units != None && units != None)
		log_fatal("Adding timestreams of unequal units");
	for (size_t i = 0; i < size(); i++)
		(*this)[i] += r[i];

	return *this;
}

G3Timestream &G3Timestream::operator -=(const G3Timestream &r)
{
	if (r.size() != size())
		log_fatal("Subtracting timestreams of unequal length");
	if (r.units != units && r.units != None && units != None)
		log_fatal("Subtracting timestreams of unequal units");
	for (size_t i = 0; i < size(); i++)
		(*this)[i] -= r[i];

	return *this;
}

G3Timestream G3Timestream::operator -(const G3Timestream &r) const
{
	G3Timestream ret(*this);

	if (r.size() != size())
		log_fatal("Subtracting timestreams of unequal length");
	if (r.units != units && r.units != None && units != None)
		log_fatal("Subtracting timestreams of unequal units");
	for (size_t i = 0; i < size(); i++)
		ret[i] = (*this)[i] - r[i];

	return ret;
}

G3Timestream G3Timestream::operator *(const G3Timestream &r) const
{
	G3Timestream ret(*this);

	if (r.size() != size())
		log_fatal("Multiplying timestreams of unequal length");
	if (r.units != units && r.units != None && units != None)
		log_fatal("Multiplying timestreams of unequal units");
	for (size_t i = 0; i < size(); i++)
		ret[i] = (*this)[i] * r[i];
	if (r.units != units)
		ret.units = r.units == None ? units : r.units;

	return ret;
}

G3Timestream G3Timestream::operator /(const G3Timestream &r) const
{
	G3Timestream ret(*this);

	if (r.size() != size())
		log_fatal("Dividing timestreams of unequal length");
	if (r.units != units && r.units != None && units != None)
		log_fatal("Dividing timestreams of unequal units");
	for (size_t i = 0; i < size(); i++)
		ret[i] = (*this)[i] / r[i];
	if (r.units == units)
		ret.units = None;

	return ret;
}

G3Timestream G3Timestream::operator +(double r) const
{
	G3Timestream ret(*this);
	for (size_t i = 0; i < size(); i++)
		ret[i] = (*this)[i] + r;
	return ret;
}

G3Timestream operator +(double l, const G3Timestream &r)
{
	return r+l;
}

G3Timestream G3Timestream::operator -(double r) const
{
	G3Timestream ret(*this);
	for (size_t i = 0; i < size(); i++)
		ret[i] = (*this)[i] - r;
	return ret;
}

G3Timestream operator -(double l, const G3Timestream &r)
{
	G3Timestream ret(r);
	for (size_t i = 0; i < r.size(); i++)
		ret[i] = l - r[i];
	return ret;
}

G3Timestream G3Timestream::operator *(double r) const
{
	G3Timestream ret(*this);
	for (size_t i = 0; i < size(); i++)
		ret[i] *= r;
	return ret;
}

G3Timestream &G3Timestream::operator *=(double r)
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] *= r;
	return *this;
}

G3Timestream &G3Timestream::operator /=(double r)
{
	for (size_t i = 0; i < size(); i++)
		(*this)[i] /= r;
	return *this;
}

G3Timestream operator *(double l, const G3Timestream &r)
{
	return r*l;
}

G3Timestream G3Timestream::operator /(double r) const
{
	G3Timestream ret(*this);
	for (size_t i = 0; i < size(); i++)
		ret[i] = (*this)[i] / r;
	return ret;
}

G3Timestream operator /(double l, const G3Timestream &r)
{
	G3Timestream ret(r);
	for (size_t i = 0; i < r.size(); i++)
		ret[i] = l/r[i];
	return ret;
}

std::string G3Timestream::Description() const
{
	std::ostringstream desc;
	desc.precision(1);
	desc << std::fixed;
	desc << size() << " samples at " << GetSampleRate()/G3Units::Hz <<
	    " Hz";
	switch (units) {
	case Counts:
		desc << " (Counts)";
		break;
	case Current:
		desc << " (Current)";
		break;
	case Power:
		desc << " (Power)";
		break;
	case Resistance:
		desc << " (Resistance)";
		break;
	case Tcmb:
		desc << " (Tcmb)";
		break;
	case Angle:
		desc << " (Angle)";
		break;
	case Distance:
		desc << " (Distance)";
		break;
	case Voltage:
		desc << " (Voltage)";
		break;
	case Pressure:
		desc << " (Pressure)";
		break;
	case FluxDensity:
		desc << " (FluxDensity)";
		break;
	case Trj:
		desc << " (Trj)";
		break;
	case Frequency:
		desc << " (Frequency)";
		break;
	default:
		break;
	}
	return desc.str();
}

template <class A> void G3TimestreamMap::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	if (v < 3) {
		std::map<std::string, G3Timestream> oldmap;
		ar & cereal::make_nvp("map", oldmap);
		for (auto i = oldmap.begin(); i != oldmap.end(); i++)
			this->insert(std::pair<std::string, G3TimestreamPtr>(
			    i->first, G3TimestreamPtr(new G3Timestream(
			    i->second))));
	} else {
		ar & cereal::make_nvp("map",
		    cereal::base_class<OrderedMap<std::string,
		    G3TimestreamPtr> >(this));
	}
	if (v < 2) {
		// Load old timestreams with start/stop in the map instead of
		// the individual timestreams.
		G3Time start, stop;
		ar & cereal::make_nvp("start", start);
		ar & cereal::make_nvp("stop", stop);
		for (auto i = begin(); i != end(); i++) {
			i->second->start = start;
			i->second->stop = stop;
		}
	}
}

std::string G3TimestreamMap::Description() const
{
	std::ostringstream desc;
	desc << "Timestreams";
	if (begin() != end())
		desc << " of " << *(begin()->second);
	desc << " from " << size() << " detectors";
	return desc.str();
}

bool G3TimestreamMap::CheckAlignment() const
{
	if (begin() == end())
		return true;

	G3Time start(begin()->second->start), stop(begin()->second->stop);
	size_t nsamps(begin()->second->size());
	for (auto i = begin(); i != end(); i++) {
		if (i->second->start.time != start.time)
			return false;
		if (i->second->stop.time != stop.time)
			return false;
		if (i->second->size() != nsamps)
			return false;
	}

	return true;
}

G3Time G3TimestreamMap::GetStartTime() const
{
	if (begin() == end())
		return G3Time();

	return begin()->second->start;
}

void G3TimestreamMap::SetStartTime(G3Time start)
{
	for (auto& ts : *this)
		ts.second->start = start;
}

G3Time G3TimestreamMap::GetStopTime() const
{
	if (begin() == end())
		return G3Time();

	return begin()->second->stop;
}

void G3TimestreamMap::SetStopTime(G3Time stop)
{
	for (auto& ts : *this)
		ts.second->stop = stop;
}

double G3TimestreamMap::GetSampleRate() const
{
	if (begin() == end())
		return 0;

	return begin()->second->GetSampleRate();
}

size_t G3TimestreamMap::NSamples() const
{
	if (begin() == end())
		return 0;

	return begin()->second->size();
}

G3Timestream::TimestreamUnits G3TimestreamMap::GetUnits() const
{
	if (begin() == end())
		return G3Timestream::None;

	return begin()->second->units;
}

void G3TimestreamMap::SetUnits(G3Timestream::TimestreamUnits units)
{
	for (auto& ts : *this)
		ts.second->units = units;
}

uint8_t G3TimestreamMap::GetFLACCompression() const
{
	if (begin() == end())
		return 0;

	return begin()->second->use_flac_;
}

uint8_t G3TimestreamMap::GetFLACBitDepth() const
{
	if (begin() == end())
		return 0;

	return begin()->second->flac_depth_;
}

void G3TimestreamMap::SetFLACCompression(int compression_level)
{
	// Check for errors
	begin()->second->SetFLACCompression(compression_level);

	for (auto& ts : *this)
		ts.second->use_flac_ = compression_level;
}

void G3TimestreamMap::SetFLACBitDepth(int bit_depth)
{
	// Check for errors
	begin()->second->SetFLACBitDepth(bit_depth);

	for (auto& ts : *this)
		ts.second->flac_depth_ = bit_depth;
}

void G3TimestreamMap::Compactify()
{
	// Check if already compacted
	bool is_compact = true;
	std::shared_ptr<void> first_root;
	void *expected_base = NULL;
	for (auto & i : *this) {
		// Not compact if any using internal storage
		if (!i.second->root_data_ref_) {
			is_compact = false;
			break;
		}
		if (!first_root)
			first_root = i.second->root_data_ref_;
		// or if they are using *different* internal storage
		if (i.second->root_data_ref_ != first_root) {
			is_compact = false;
			break;
		}
		// or if they are not in order
		size_t element_size = 0;
		switch (i.second->data_type_) {
		case G3Timestream::TS_DOUBLE:
		case G3Timestream::TS_INT64:
			element_size = 8;
			break;
		case G3Timestream::TS_FLOAT:
		case G3Timestream::TS_INT32:
			element_size = 4;
			break;
		default:
			log_fatal("Unknown timestream datatype %d", i.second->data_type_);
		}
		if (expected_base != NULL && i.second->data_ != expected_base) {
			is_compact = false;
			break;
		}
		expected_base = (uint8_t *)i.second->data_ +
		    element_size*i.second->size();
	}
	if (is_compact || size() == 0)
		return;

	// Check if timestreams aligned; if not, they can't be compacted
	if (!CheckAlignment())
		throw std::runtime_error("Cannot compactify a timestream map "
		    "with non-aligned timestreams.");

	// Check if all timestreams have the same data type
	for (auto & i : *this) {
		if (i.second->data_type_ != begin()->second->data_type_)
			throw std::runtime_error("Cannot compactify a "
			    "timestream map with mixed data types.");
	}

	// Now let's do the compactification.
	std::shared_ptr<void> root;
	void *base_ptr;
	size_t row_bytes;
	switch (begin()->second->data_type_) {
	case G3Timestream::TS_DOUBLE: {
		std::vector<double> *data = new std::vector<double>(size() *
		    begin()->second->size());
		root = std::shared_ptr<std::vector<double> >(data);
		base_ptr = &(*data)[0];
		row_bytes = begin()->second->size() * sizeof((*data)[0]);
		break;
	}
	case G3Timestream::TS_FLOAT: {
		std::vector<float> *data = new std::vector<float>(size() *
		    begin()->second->size());
		root = std::shared_ptr<std::vector<float> >(data);
		base_ptr = &(*data)[0];
		row_bytes = begin()->second->size() * sizeof((*data)[0]);
		break;
	}
	case G3Timestream::TS_INT32: {
		std::vector<int32_t> *data = new std::vector<int32_t>(size() *
		    begin()->second->size());
		root = std::shared_ptr<std::vector<int32_t> >(data);
		base_ptr = &(*data)[0];
		row_bytes = begin()->second->size() * sizeof((*data)[0]);
		break;
	}
	case G3Timestream::TS_INT64: {
		std::vector<int64_t> *data = new std::vector<int64_t>(size() *
		    begin()->second->size());
		root = std::shared_ptr<std::vector<int64_t> >(data);
		base_ptr = &(*data)[0];
		row_bytes = begin()->second->size() * sizeof((*data)[0]);
		break;
	}
	default:
		log_fatal("Unknown timestream datatype %d", begin()->second->data_type_);
	}

	for (auto & i : *this) {
		memcpy(base_ptr, i.second->data_, row_bytes);
		if (i.second->buffer_)
			delete i.second->buffer_;
		i.second->buffer_ = NULL;
		i.second->root_data_ref_ = root;
		i.second->data_ = base_ptr;
		// len unchanged

		base_ptr = (uint8_t *)base_ptr + row_bytes;
	}
}

G3_SPLIT_SERIALIZABLE_CODE(G3Timestream);
G3_SERIALIZABLE_CODE(G3TimestreamMap);

static void
G3Timestream_assert_congruence(const G3Timestream &a, const G3Timestream &b)
{

	if (b.size() != a.size())
		log_fatal("Timestreams of unequal length");
	if (a.units != b.units && a.units != G3Timestream::None &&
	    b.units != G3Timestream::None)
		log_fatal("Timestreams of unequal units");
	if (a.start.time != b.start.time)
		log_fatal("Timestreams start at different times");
	if (a.stop.time != b.stop.time)
		log_fatal("Timestreams stop at different times");
}

static G3TimestreamPtr
G3Timestream_getslice(const G3Timestream &a, py::slice slice)
{
	int start(0), stop(a.size()), step(1);
	double sample_spacing = 1./a.GetSampleRate();

	// Normalize and check slice boundaries
	if (slice.start().ptr() != Py_None)
		start = py::extract<int>(slice.start())();
	if (slice.stop().ptr() != Py_None)
		stop = py::extract<int>(slice.stop())();
	if (slice.step().ptr() != Py_None)
		step = py::extract<int>(slice.step())();

	if (start < 0)
		start = a.size() + start;
	if (stop < 0)
		stop = a.size() + stop;

	if (stop > (int)a.size())
		stop = a.size();
	if (step > (int)a.size())
		step = a.size();

	if (start >= (int)a.size() || start < 0)
		log_fatal("Start index %d out of range", start);
	if (stop < 0)
		log_fatal("Stop index %d out of range", stop);
	if (step <= 0)
		log_fatal("Step index %d out of range", step);
	if (start >= stop)
		log_fatal("Start index %d >= stop index %d", start, stop);

	// Get stop index corresponding to step parameter
	stop = start + ((stop - start + (step - 1))/step)*step;

	// Build new TS
	G3TimestreamPtr out(new G3Timestream((stop - start)/step));
	out->units = a.units;
	out->start.time = a.start.time + G3TimeStamp(start*sample_spacing);
	out->stop.time = a.start.time +
	    G3TimeStamp((stop - step)*sample_spacing);

	for (int i = start, j = 0; i < stop; i += step, j++)
		(*out)[j] = a[i];

	return out;
}

static auto
G3Timestream_elapsed(const G3Timestream &a)
{
	if (!a.size())
		return G3VectorIntPtr();

	double sample_spacing = 1. / a.GetSampleRate();

	G3VectorIntPtr v(new G3VectorInt);
	v->reserve(a.size());

	for (size_t i = 0; i < a.size(); i++)
		v->push_back(i * sample_spacing);

	return v;
}

static auto
G3TimestreamMap_elapsed(const G3TimestreamMap &a)
{
	if (!a.size())
		return G3VectorIntPtr();

	return G3Timestream_elapsed(*a.begin()->second);
}

static auto
G3Timestream_times(const G3Timestream &a)
{
	if (!a.size())
		return G3VectorTimePtr();

	double sample_spacing = 1. / a.GetSampleRate();

	auto v = G3VectorTimePtr(new G3VectorTime);
	v->reserve(a.size());

	for (size_t i = 0; i < a.size(); i++)
		v->push_back(G3Time(a.start.time + (uint64_t)(i * sample_spacing)));

	return v;
}

static auto
G3TimestreamMap_times(const G3TimestreamMap &a)
{
	if (!a.size())
		return G3VectorTimePtr();

	return G3Timestream_times(*a.begin()->second);
}


class G3Timestream::G3TimestreamPythonHelpers
{
public:
SET_LOGGER("G3Timestream");

static DataType
get_ts_dtype(const Py_buffer &info)
{
	std::string format = check_buffer_format(info.format);

	if (format == "d") {
		return TS_DOUBLE;
	} else if (format == "f") {
		return TS_FLOAT;
#ifdef __LP64__
	} else if (format == "i") {
#else
	} else if (format == "i" || format == "l") {
#endif
		assert(info.itemsize == sizeof(int32_t));
		return TS_INT32;
#ifdef __LP64__
	} else if (format == "q" || format == "l") {
#else
	} else if (format == "q") {
#endif
		assert(info.itemsize == sizeof(int64_t));
		return TS_INT64;
	} else {
		throw py::type_error(std::string("Unsupported data type: ") + info.format);
	}
}

static int
G3TimestreamMap_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;
	view->internal = NULL;
	view->suboffsets = NULL;
	view->buf = NULL;

	py::handle<> self(py::borrowed(obj));
	py::object selfobj(self);
	py::extract<G3TimestreamMapPtr> ext(selfobj);
	if (!ext.check()) {
		PyErr_SetString(PyExc_ValueError, "Invalid timestream");
		view->obj = NULL;
		return -1;
	}
	G3TimestreamMapPtr ts = ext();
	if (!ts->CheckAlignment()) {
		PyErr_SetString(PyExc_BufferError, "Timestream map is not "
		    "aligned, cannot cast to a 2D array.");
		view->obj = NULL;
		return -1;
	}
	if (ts->size() == 0) {
		PyErr_SetString(PyExc_BufferError, "Timestream map is empty.");
		view->obj = NULL;
		return -1;
	}

	if ((flags & PyBUF_F_CONTIGUOUS) == PyBUF_F_CONTIGUOUS) {
		PyErr_SetString(PyExc_BufferError, "Cannot provide FORTRAN "
		    "contiguous buffer.");
		view->obj = NULL;
		return -1;
	}

	// Get everything into a 2D array
	ts->Compactify();

	view->obj = obj;
	view->readonly = 0;
	switch (ts->begin()->second->data_type_) {
	case TS_DOUBLE:
		view->format = (char *)"d";
		view->itemsize = sizeof(double);
		break;
	case TS_FLOAT:
		view->format = (char *)"f";
		view->itemsize = sizeof(float);
		break;
	case TS_INT32:
		view->format = (char *)"i";
		view->itemsize = sizeof(int32_t);
		break;
	case TS_INT64:
		view->format = (char *)"q";
		view->itemsize = sizeof(int64_t);
		break;
	default:
		PyErr_SetString(PyExc_BufferError, "Unsupported data type.");
		view->obj = NULL;
		return -1;
	}
	if (!(flags & PyBUF_FORMAT))
		view->format = NULL;
	view->ndim = 2;
	view->len = ts->size() * ts->begin()->second->size() * view->itemsize;

	view->shape = new Py_ssize_t[2];
	view->shape[0] = ts->size();
	view->shape[1] = ts->begin()->second->size();

	view->buf = ts->begin()->second->data_;
	view->readonly = 0;

	view->strides = new Py_ssize_t[2];
	view->strides[0] = ts->begin()->second->size()*view->itemsize;
	view->strides[1] = view->itemsize; 

	Py_INCREF(obj);

	return 0;
}

static void
G3TimestreamMap_relbuffer(PyObject *obj, Py_buffer *view)
{
	if (view->strides != NULL)
		delete [] view->strides;
	if (view->shape != NULL)
		delete [] view->shape;
}

static G3TimestreamPtr
timestream_from_iterable(py::object v,
    G3Timestream::TimestreamUnits units = G3Timestream::None)
{
	// Sometimes the explicit copy constructor does not get
	// priority from python, so do what that should have done.
	py::extract<G3TimestreamConstPtr> was_ts_already(v);
	if (was_ts_already.check())
		return G3TimestreamPtr(new G3Timestream(*was_ts_already()));

	// That out of the way, move on to the generic numpy-ish case
	G3TimestreamPtr x;
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_CONTIG_RO) != -1) {
		if (strcmp(view.format, "d") == 0) {
			x = G3TimestreamPtr(new G3Timestream((double *)view.buf,
			    (double *)view.buf + view.len/sizeof(double)));
		} else if (strcmp(view.format, "f") == 0) {
			x = G3TimestreamPtr(new G3Timestream());
			delete x->buffer_; x->buffer_ = NULL;
			x->data_type_ = G3Timestream::TS_FLOAT;
			float *data = new float[view.len/sizeof(float)];
			x->root_data_ref_ = std::shared_ptr<float[]>(data);
			x->data_ = data;
			x->len_ = view.len/sizeof(float);
			memcpy(data, view.buf, view.len);
#ifdef __LP64__
		} else if (strcmp(view.format, "i") == 0) {
#else
		} else if (strcmp(view.format, "i") == 0 || strcmp(view.format, "l") == 0) {
#endif
			x = G3TimestreamPtr(new G3Timestream());
			delete x->buffer_; x->buffer_ = NULL;
			x->data_type_ = G3Timestream::TS_INT32;
			int32_t *data = new int32_t[view.len/sizeof(int32_t)];
			x->root_data_ref_ = std::shared_ptr<int32_t[]>(data);
			x->data_ = data;
			x->len_ = view.len/sizeof(int32_t);
			memcpy(data, view.buf, view.len);
			assert(view.itemsize == sizeof(int32_t));
#ifdef __LP64__
		} else if (strcmp(view.format, "q") == 0 || strcmp(view.format, "l") == 0) {
#else
		} else if (strcmp(view.format, "q") == 0) {
#endif
			x = G3TimestreamPtr(new G3Timestream());
			delete x->buffer_; x->buffer_ = NULL;
			x->data_type_ = G3Timestream::TS_INT64;
			int64_t *data = new int64_t[view.len/sizeof(int64_t)];
			x->root_data_ref_ = std::shared_ptr<int64_t[]>(data);
			x->data_ = data;
			x->len_ = view.len/sizeof(int64_t);
			memcpy(data, view.buf, view.len);
			assert(view.itemsize == sizeof(int64_t));
		} else {
			// We could add more types, but why do that?
			// Let Python do the work for obscure cases
			std::vector<double> data;
			py::container_utils::extend_container(data, v);
			x = G3TimestreamPtr(new G3Timestream(data.begin(), data.end()));
		}
		PyBuffer_Release(&view);
	} else {
		PyErr_Clear();
		std::vector<double> data;
		py::container_utils::extend_container(data, v);
		x = G3TimestreamPtr(new G3Timestream(data.begin(), data.end()));
	}

	x->units = units;

	return x;
}

struct PyBufferOwner {
	PyBufferOwner(Py_buffer view) : v(view) {}
	~PyBufferOwner() { PyBuffer_Release(&v); }
	Py_buffer v;
};

static G3TimestreamMapPtr
G3TimestreamMap_from_numpy(std::vector<std::string> keys,
    py::object data, G3Time start, G3Time stop,
    G3Timestream::TimestreamUnits units, int compression_level,
    bool copy_data, int bit_depth)
{
	G3TimestreamMapPtr x(new G3TimestreamMap);

	std::shared_ptr<PyBufferOwner> v;

	{
	Py_buffer view;
	if (PyObject_GetBuffer(data.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) == -1) {
		PyErr_Clear();
		
		// Maybe it's an indexable iterable? See if processing it
		// that way throws an exception.
		for (size_t i = 0; i < keys.size(); i++) {
			G3TimestreamPtr next =
			    timestream_from_iterable(data[i], units);
			next->start = start;
			next->stop = stop;
			next->SetFLACCompression(compression_level);
			next->SetFLACBitDepth(bit_depth);
			(*x)[keys[i]] = next;
		}
		return x;
	}

	// Set up an auto-deleter for view when the function exits now
	// that it has been allocated successfully. This is a shared pointer
	// so that it ownership can be moved to the timestreams if needed.
	v = std::make_shared<PyBufferOwner>(view);
	}

	if (keys.size() != (size_t)v->v.shape[0])
		throw py::index_error("Number of keys does not "
		    "match number of rows in data structure.");

	if (v->v.ndim != 2)
		throw py::value_error("Array must be two-dimensional.");

	G3Timestream templ;
	templ.units = units;
	templ.start = start;
	templ.stop = stop;
	templ.SetFLACCompression(compression_level);
	templ.SetFLACBitDepth(bit_depth);
	templ.data_type_ = get_ts_dtype(v->v);

	uint8_t *buf;
	std::shared_ptr<void> data_ref;
	size_t len = v->v.shape[1];
	ptrdiff_t step = v->v.strides[0];
	if (copy_data) {
		buf = new uint8_t[v->v.len];
		data_ref = std::shared_ptr<uint8_t[]>(buf);
		memcpy(buf, v->v.buf, v->v.len);
		v.reset(); // Release Python Buffer view
	} else {
		buf = (uint8_t *)v->v.buf;
		data_ref = v; // Keep view around
	}

	for (auto &i : keys) {
		G3TimestreamPtr next(new G3Timestream(templ));
		delete next->buffer_;
		next->buffer_ = NULL;
		next->len_ = len;
		next->data_ = buf;
		next->root_data_ref_ = data_ref;
		(*x)[i] = next;

		buf += step;
	}

	return x;
}

static int
G3Timestream_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	py::handle<> self(py::borrowed(obj));
	py::object selfobj(self);
	py::extract<G3TimestreamPtr> ext(selfobj);
	if (!ext.check()) {
		PyErr_SetString(PyExc_ValueError, "Invalid timestream");
		view->obj = NULL;
		return -1;
	}
	G3TimestreamPtr ts = ext();
	view->obj = obj;
	view->buf = ts->data_;
	view->readonly = 0;
	switch (ts->data_type_) {
	case TS_DOUBLE:
		view->format = (char *)"d";
		view->itemsize = sizeof(double);
		break;
	case TS_FLOAT:
		view->format = (char *)"f";
		view->itemsize = sizeof(float);
		break;
	case TS_INT32:
		view->format = (char *)"i";
		view->itemsize = sizeof(int32_t);
		break;
	case TS_INT64:
		view->format = (char *)"q";
		view->itemsize = sizeof(int64_t);
		break;
	default:
		PyErr_SetString(PyExc_ValueError, "Unsupported data type.");
		view->obj = NULL;
		return -1;
	}
	if (!(flags & PyBUF_FORMAT))
		view->format = NULL;
	view->len = ts->len_ * view->itemsize;
	view->ndim = 1;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION >= 3)
	// Abuse internal pointer in the absence of smalltable. This is safe
	// on all architectures except MIPS N32.
	view->internal = (void *)ts->size();
	view->shape = (Py_ssize_t *)(&view->internal);
#else
	view->smalltable[0] = ts->size();
	view->shape = &view->smalltable[0];
	view->internal = NULL;
#endif
	view->strides = &view->itemsize;
	view->suboffsets = NULL;

	// Try to hold onto our collective hats. This is still very dangerous if
	// the G3Timestream's underlying vector is resized.
	Py_INCREF(obj);

	return 0;
}

};

static py::tuple
G3Timestream_shape(const G3Timestream &r)
{
	return py::make_tuple(r.size());
}

static int
G3Timestream_ndim(const G3Timestream &r)
{
	return 1;
}

static PyBufferProcs timestream_bufferprocs;
static PyBufferProcs timestreammap_bufferprocs;

PYBINDINGS("core", scope) {
	register_enum<G3Timestream::TimestreamUnits, G3Timestream::None>(scope,
	  "G3TimestreamUnits",
	  "Unit scheme for timestreams and maps. Designates different classes "
	  "of units (power, current, on-sky temperature) rather than choices "
	  "of unit within a class (watts vs. horsepower, or K vs. uK), "
	  "transformations between which are handled by core.G3Units.")
	    .value("None",  G3Timestream::None)
	    .value("Counts",  G3Timestream::Counts)
	    .value("Current",  G3Timestream::Current)
	    .value("Power",  G3Timestream::Power)
	    .value("Resistance",  G3Timestream::Resistance)
	    .value("Tcmb",  G3Timestream::Tcmb)
	    .value("Angle",  G3Timestream::Angle)
	    .value("Distance",  G3Timestream::Distance)
	    .value("Voltage",  G3Timestream::Voltage)
	    .value("Pressure",  G3Timestream::Pressure)
	    .value("FluxDensity",  G3Timestream::FluxDensity)
	    .value("Trj",  G3Timestream::Trj)
	    .value("Frequency",  G3Timestream::Frequency)
	;

	py::object ts =
	register_frameobject<G3Timestream>(scope, "G3Timestream",
	   "Detector timestream. "
	   "Includes a units field and start and stop times. Can otherwise be "
	   "treated as a numpy array with a float64 dtype. Conversions to and "
           "from such arrays (e.g. with numpy.asarray) are fast. Note that a "
           "numpy array constructed from a timestream will share a memory "
           "buffer: changes to the array affect the timestream and vice versa. "
	   "Most binary timestream arithmetic operations (+, -) check that the "
	   "units and start/stop times are congruent.")
	    .def(py::init<>())
	    .def("__init__", py::make_constructor(G3Timestream::G3TimestreamPythonHelpers::timestream_from_iterable, py::default_call_policies(), (py::arg("data"), py::arg("units") = G3Timestream::TimestreamUnits::None)), "Create a timestream from a numpy array or other numeric python iterable")
	    .def("SetFLACCompression", &G3Timestream::SetFLACCompression,
	      "Pass True to turn on FLAC compression when serialized. "
	      "FLAC compression only works if the timestream is in units of "
	      "counts.")
	    .def("SetFLACBitDepth", &G3Timestream::SetFLACBitDepth,
	      "Change the bit depth for FLAC compression, may be 24 or 32 "
	      "(default, requires version 1.4+).")
	    .def_readwrite("units", &G3Timestream::units,
	      "Units of the data in the timestream, stored as one of the "
	      "members of core.G3TimestreamUnits.")
	    .def_readwrite("start", &G3Timestream::start,
	      "Time of the first sample in the time stream")
	    .def_readwrite("stop", &G3Timestream::stop,
	      "Time of the final sample in the timestream")
	    .add_property("sample_rate", &G3Timestream::GetSampleRate,
	      "Computed sample rate of the timestream.")
	    .add_property("n_samples", &G3Timestream::size,
	      "Number of samples in the timestream. Equivalent to len(ts)")
	    .add_property("compression_level", &G3Timestream::GetFLACCompression,
	      &G3Timestream::SetFLACCompression, "Level of FLAC compression used for this timestream. "
	      "This can only be non-zero if the timestream is in units of counts.")
	    .add_property("bit_depth", &G3Timestream::GetFLACBitDepth,
	      &G3Timestream::SetFLACBitDepth, "Bit depth of FLAC compression used for this timestream.")
	    .def("_assert_congruence", &G3Timestream_assert_congruence,
	      "log_fatal() if units, length, start, or stop times do not match")
	    .def("_cxxslice", &G3Timestream_getslice, "Slice-only __getitem__")
	    .add_property("elapsed", &G3Timestream_elapsed,
	      "Compute elapsed time array for samples")
	    .add_property("times", &G3Timestream_times,
	      "Compute time vector for samples")
	    .def("__len__", &G3Timestream::size)
	    .add_property("shape", &G3Timestream_shape)
	    .add_property("ndim", &G3Timestream_ndim)
	    // Operators bound in python through numpy
	;

	// Add buffer protocol interface
	PyTypeObject *tsclass = (PyTypeObject *)ts.ptr();
	timestream_bufferprocs.bf_getbuffer = G3Timestream::G3TimestreamPythonHelpers::G3Timestream_getbuffer;
	tsclass->tp_as_buffer = &timestream_bufferprocs;
#if PY_MAJOR_VERSION < 3
	tsclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	py::object tsm =
	register_frameobject<G3TimestreamMap>(scope, "G3TimestreamMap",
	    "Collection of timestreams indexed by logical detector ID")
	    .def(py::init<>())
	    .def("__init__", py::make_constructor(G3Timestream::G3TimestreamPythonHelpers::G3TimestreamMap_from_numpy, 
	         py::default_call_policies(),
	         (py::arg("keys"), py::arg("data"), py::arg("start")=G3Time(0),
	          py::arg("stop")=G3Time(0), py::arg("units") = G3Timestream::TimestreamUnits::None,
	          py::arg("compression_level") = 0, py::arg("copy_data") = true, py::arg("bit_depth") = 32)),
	         "Create a timestream map from a numpy array or other numeric python iterable. "
	         "Each row of the 2D input array will correspond to a single timestream, with "
	         "the key set to the correspondingly-indexed entry of <keys>. If <copy_data> "
	         "is True (default), the data will be copied into the output data structure. "
	         "If False, the timestream map will provide a view into the given numpy array.")
	    .def(py::std_map_indexing_suite<G3TimestreamMap, true>())
	    .def("CheckAlignment", &G3TimestreamMap::CheckAlignment)
	    .def("Compactify", &G3TimestreamMap::Compactify,
	       "If member timestreams are stored non-contiguously, repack all "
	       "data into a contiguous block. Requires timestreams be aligned "
	       "and the same data type. Done implicitly by numpy.asarray().")
	    .def("SetFLACCompression", &G3TimestreamMap::SetFLACCompression,
	      "Pass True to turn on FLAC compression when serialized. "
	      "FLAC compression only works if the timestreams are in units of "
	      "counts.")
	    .def("SetFLACBitDepth", &G3TimestreamMap::SetFLACBitDepth,
	      "Change the bit depth for FLAC compression, may be 24 or 32 "
	      "(default, requires version 1.4+).")
	    .add_property("start", &G3TimestreamMap::GetStartTime,
	      &G3TimestreamMap::SetStartTime,
	      "Time of the first sample in the time stream")
	    .add_property("stop", &G3TimestreamMap::GetStopTime,
	      &G3TimestreamMap::SetStopTime,
	      "Time of the final sample in the time stream")
	    .add_property("sample_rate", &G3TimestreamMap::GetSampleRate,
	      "Computed sample rate of the timestream.")
	    .add_property("n_samples", &G3TimestreamMap::NSamples,
	      "Number of samples in the timestream. Equivalent to the length "
	      "of one of the timestreams.")
	    .add_property("units", &G3TimestreamMap::GetUnits,
	      &G3TimestreamMap::SetUnits,
	      "Units of the data in the timestream, stored as one of the "
	      "members of core.G3TimestreamUnits.")
	    .add_property("compression_level", &G3TimestreamMap::GetFLACCompression,
	      &G3TimestreamMap::SetFLACCompression,
	      "Level of FLAC compression used for this timestream map. "
	      "This can only be non-zero if the timestream is in units of counts.")
	    .add_property("bit_depth", &G3TimestreamMap::GetFLACBitDepth,
	      &G3TimestreamMap::SetFLACBitDepth,
	      "Bit depth of FLAC compression used for this timestream map.")
	    .add_property("elapsed", &G3TimestreamMap_elapsed,
	      "Compute elapsed time array for samples")
	    .add_property("times", &G3TimestreamMap_times,
	      "Compute time vector for samples")
	;

	// Add buffer protocol interface
	PyTypeObject *tsmclass = (PyTypeObject *)tsm.ptr();
	timestreammap_bufferprocs.bf_getbuffer =
	    G3Timestream::G3TimestreamPythonHelpers::G3TimestreamMap_getbuffer;
	timestreammap_bufferprocs.bf_releasebuffer =
	    G3Timestream::G3TimestreamPythonHelpers::G3TimestreamMap_relbuffer;
	tsmclass->tp_as_buffer = &timestreammap_bufferprocs;
#if PY_MAJOR_VERSION < 3
	tsmclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif
}

