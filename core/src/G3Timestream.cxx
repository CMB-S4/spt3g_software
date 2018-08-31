#include <pybindings.h>
#include <serialization.h>
#include <G3Timestream.h>
#include <std_map_indexing_suite.hpp>
#include <boost/python/slice.hpp>

#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>

#ifdef G3_HAS_FLAC
#include <FLAC/stream_encoder.h>
#include <cmath>

template<typename A>
struct FlacDecoderCallbackArgs {
	A *inbuf;
	std::vector<double> *outbuf;
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

		if (units != Counts)
			log_fatal("Cannot use FLAC on non-counts timestreams");

		// Copy to 24-bit integers
		inbuf.resize(size());
		for (size_t i = 0; i < size(); i++)
			inbuf[i] = ((int32_t((*this)[i]) & 0x00ffffff) << 8)
			    >> 8;
		chanmap[0] = &inbuf[0];

		// Mark bad samples using an out-of-band signal. Since we
		// convert to 24-bit integers going into FLAC, which have no
		// out-of-range values for signalling, this requires special
		// care. Usually a timestream is either all-valid or
		// all-invalid, so signal that with a single byte flag. In the
		// rare case that only some samples are valid, store a
		// validity mask.
		std::vector<bool> nanbuf(size(), false);
		for (size_t i = 0; i < size(); i++) {
			if (!std::isfinite((*this)[i])) {
				nans++;
				nanbuf[i] = true;
				inbuf[i] = 0;
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
		// XXX: should assert if high-order 8 bits are not clear
		FLAC__stream_encoder_set_bits_per_sample(encoder, 24);
		FLAC__stream_encoder_set_compression_level(encoder, use_flac_);
		FLAC__stream_encoder_init_stream(encoder,
		    flac_encoder_write_cb, NULL, NULL, NULL, (void*)(&outbuf));
		FLAC__stream_encoder_process (encoder, chanmap, inbuf.size());
		FLAC__stream_encoder_finish(encoder);
		FLAC__stream_encoder_delete(encoder);

		ar & cereal::make_nvp("data", outbuf);
	} else {
#endif
		ar & cereal::make_nvp("data",
		    cereal::base_class<std::vector<double> >(this));
#ifdef G3_HAS_FLAC
	}
#endif
}

template <class A> void G3Timestream::load(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("units", units);
	if (v > 1) {
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
		callback.outbuf = this;
		callback.pos = 0;

		if (units != Counts)
			log_fatal("Cannot use FLAC on non-counts timestreams");

		ar & cereal::make_nvp("nanflag", nanflag);
		if (nanflag == SomeNan)
			ar & cereal::make_nvp("nanmask", nanbuf);

		ar & cereal::make_size_tag(callback.nbytes);

		// Typical compression ratio: N bytes in input = N samples
		reserve(callback.nbytes);

		FLAC__StreamDecoder *decoder = FLAC__stream_decoder_new();
		FLAC__stream_decoder_init_stream(decoder,
		    flac_decoder_read_cb<A>, NULL, NULL, NULL, NULL,
		    flac_decoder_write_cb<A>, NULL, flac_decoder_error_cb,
		    (void*)(&callback));
		FLAC__stream_decoder_process_until_end_of_stream(decoder);
		FLAC__stream_decoder_finish(decoder);
		FLAC__stream_decoder_delete(decoder);

		// Apply NaN mask
		if (nanflag == AllNan) {
			for (size_t i = 0; i < size(); i++)
				(*this)[i] = NAN;
		} else if (nanflag == SomeNan) {
			for (size_t i = 0; i < size(); i++)
				if (nanbuf[i])
					(*this)[i] = NAN;
		}

#else
		log_fatal("Trying to read FLAC-compressed timestreams but built without FLAC support");
#endif
	} else {
		ar & cereal::make_nvp("data",
		    cereal::base_class<std::vector<double> >(this));
	}
}

double G3Timestream::GetSampleRate() const
{
	return double(size() - 1)/double(stop.time - start.time);
}

void G3Timestream::SetFLACCompression(int use_flac)
{

#ifdef G3_HAS_FLAC
	if (use_flac != 0 && units != Counts)
		log_fatal("Cannot use FLAC on non-counts timestreams");

	use_flac_ = use_flac;
#else
	if (use_flac != 0)
		log_fatal("Built without FLAC support");
#endif
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
	ret.units = None;

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

static
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

static
G3Timestream operator -(double l, const G3Timestream &r)
{
	return (r-l)*(-1);
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

static
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

static
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
		    cereal::base_class<std::map<std::string,
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
	desc << "Timestreams from " << size() << " detectors";
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

static void timestream_map_set_start_time(G3TimestreamMap &map, G3Time start)
{
	for (auto i = map.begin(); i != map.end(); i++)
		i->second->start = start;
}

G3Time G3TimestreamMap::GetStopTime() const
{
	if (begin() == end())
		return G3Time();

	return begin()->second->stop;
}

static void timestream_map_set_stop_time(G3TimestreamMap &map, G3Time stop)
{
	for (auto i = map.begin(); i != map.end(); i++)
		i->second->stop = stop;
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

G3_SPLIT_SERIALIZABLE_CODE(G3Timestream);
G3_SERIALIZABLE_CODE(G3TimestreamMap);

namespace {
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

	boost::python::handle<> self(boost::python::borrowed(obj));
	boost::python::object selfobj(self);
	G3TimestreamMapPtr ts =
	    boost::python::extract<G3TimestreamMapPtr>(selfobj)();
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
#if 0
	if ((flags & (PyBUF_WRITABLE | PyBUF_ANY_CONTIGUOUS)) ==
	    (PyBUF_WRITABLE | PyBUF_ANY_CONTIGUOUS)) {
#else
	if ((flags & PyBUF_WRITABLE) == PyBUF_WRITABLE) { // See XXX note below
#endif
		PyErr_SetString(PyExc_BufferError, "Cannot provide writable "
		    "contiguous buffer.");
		view->obj = NULL;
		return -1;
	}
	if ((flags & PyBUF_F_CONTIGUOUS) == PyBUF_F_CONTIGUOUS) {
		PyErr_SetString(PyExc_BufferError, "Cannot provide FORTRAN "
		    "contiguous buffer.");
		view->obj = NULL;
		return -1;
	}

	view->obj = obj;
	view->len = ts->size() * ts->begin()->second->size() * sizeof(double);
	view->readonly = 0;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;
	view->ndim = 2;

	view->shape = new Py_ssize_t[2];
	view->shape[0] = ts->size();
	view->shape[1] = ts->begin()->second->size();

	if (false && (flags & PyBUF_INDIRECT) == PyBUF_INDIRECT) {
		// XXX: code path disabled as a result of a bug in numpy:
		// it calls PyMemoryView_FromObject, which implicitly
		// sets PyBUF_INDIRECT, but then has no code to deal
		// with suboffsets. The code in this branch does work, but
		// breaks numpy.
		
		// Provide an "indirect" buffer in which the buffer
		// is an array of pointers to buffers. The API is really
		// weird for this.
		//
		// NB: This *must* get used before any changes to the provider!
		// We do *not* increment reference counts on timestreams, nor
		// do we even have a way to freeze the buffer locations.

		view->strides = new Py_ssize_t[2];
		view->strides[0] = sizeof(double *);
		view->strides[1] = sizeof(double);

		view->suboffsets = new Py_ssize_t[2];
		view->suboffsets[0] = 0;
		view->suboffsets[1] = -1;

		view->buf = malloc(sizeof(double *)*ts->size());

		int j = 0;
		for (auto i : (*ts))
			((double **)(view->buf))[j++] = &(*i.second)[0];
	} else {
		// To honor a contiguous buffer request, make a copy of the
		// data. This violates the spirit of the buffer protocol
		// slightly, but both simplifies the API and allows a
		// potential faster memory copy than iterating over the
		// map in Python would.

		view->buf = malloc(view->len);
		view->readonly = 1;

		view->strides = new Py_ssize_t[2];
		view->strides[0] = ts->begin()->second->size()*view->itemsize;
		view->strides[1] = view->itemsize; 

		int j = 0;
		for (auto i : (*ts)) {
			memcpy((int8_t *)view->buf + j*view->strides[0],
			    &(*i.second)[0], view->strides[0]);
			j++;
		}

		view->suboffsets = NULL;
		view->internal = view->buf;
	}

	// Try to hold onto our collective hats. This is still very dangerous if
	// the G3Timestream's underlying vector is resized in the case that the
	// buffer is not a copy.
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
	if (view->suboffsets != NULL)
		delete [] view->suboffsets;
	if (view->buf != NULL)
		free(view->buf);
}
}

static PyBufferProcs timestreammap_bufferprocs;

namespace {

SET_LOGGER("G3Timestream");

G3TimestreamPtr
timestream_from_iterable(boost::python::object v,
    G3Timestream::TimestreamUnits units = G3Timestream::None)
{
	// Sometimes the explicit copy constructor does not get
	// priority from python, so do what that should have done.
	boost::python::extract<G3TimestreamConstPtr> was_ts_already(v);
	if (was_ts_already.check())
		return G3TimestreamPtr(new G3Timestream(*was_ts_already()));

	// That out of the way, move on to the generic numpy-ish case
	G3TimestreamPtr x(new G3Timestream);
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (strcmp(view.format, "d") == 0) {
			x->insert(x->begin(), (double *)view.buf,
			    (double *)view.buf + view.len/sizeof(double));
		} else if (strcmp(view.format, "f") == 0) {
			x->resize(view.len/sizeof(float));
			for (size_t i = 0; i < view.len/sizeof(float); i++)
				(*x)[i] = ((float *)view.buf)[i];
		} else {
			// We could add more types, but why do that?
			// Let Python do the work for obscure cases
			boost::python::container_utils::extend_container(*x, v);
		}
		PyBuffer_Release(&view);
	} else {
		PyErr_Clear();
		boost::python::container_utils::extend_container(*x, v);
	}

	x->units = units;

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

	boost::python::handle<> self(boost::python::borrowed(obj));
	boost::python::object selfobj(self);
	G3TimestreamPtr ts = boost::python::extract<G3TimestreamPtr>(selfobj)();
	view->obj = obj;
	view->buf = (void*)&(*ts)[0];
	view->len = ts->size() * sizeof(double);
	view->readonly = 0;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;
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

static int
G3Timestream_nsamples(const G3Timestream &r)
{
	return r.size();
}

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
G3Timestream_getslice(const G3Timestream &a, boost::python::slice slice)
{
	using namespace boost::python;
	int start(0), stop(a.size()), step(1);
	double sample_spacing = 1./a.GetSampleRate();

	// Normalize and check slice boundaries
	if (slice.start().ptr() != Py_None)
		start = extract<int>(slice.start())();
	if (slice.stop().ptr() != Py_None)
		stop = extract<int>(slice.stop())();
	if (slice.step().ptr() != Py_None)
		step = extract<int>(slice.step())();

	if (start < 0)
		start = a.size() + start;
	if (stop < 0)
		stop = a.size() + stop;

	if (start >= a.size() || start < 0)
		log_fatal("Start index %d out of range", start);
	if (stop > a.size() || stop < 0)
		log_fatal("Stop index %d out of range", stop);
	if (step >= a.size() || step <= 0)
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
}

static PyBufferProcs timestream_bufferprocs;

PYBINDINGS("core") {
	namespace bp = boost::python;

	bp::object ts =
	  EXPORT_FRAMEOBJECT(G3Timestream, init<>(), "Detector timestream. "
	   "Includes a units field and start and stop times. Can otherwise be "
	   "treated as a numpy array with a float64 dtype. Conversions to and "
           "from such arrays (e.g. with numpy.asarray) are fast. Note that a "
           "numpy array constructed from a timestream will share a memory "
           "buffer: changes to the array affect the timestream and vice versa. "
	   "Most binary timestream arithmetic operations (+, -) check that the "
	   "units and start/stop times are congruent.")
	    .def("__init__", bp::make_constructor(timestream_from_iterable, bp::default_call_policies(), (bp::arg("data"), bp::arg("units") = G3Timestream::TimestreamUnits::None)), "Create a timestream from a numpy array or other numeric python iterable")
	    .def("SetFLACCompression", &G3Timestream::SetFLACCompression,
	      "Pass True to turn on FLAC compression when serialized. "
	      "FLAC compression only works if the timestream is in units of "
	      "counts.")
	    .def_readwrite("units", &G3Timestream::units,
	      "Units of the data in the timestream, stored as one of the "
	      "members of core.G3TimestreamUnits.")
	    .def_readwrite("start", &G3Timestream::start,
	      "Time of the first sample in the time stream")
	    .def_readwrite("stop", &G3Timestream::stop,
	      "Time of the final sample in the timestream")
	    .add_property("sample_rate", &G3Timestream::GetSampleRate,
	      "Computed sample rate of the timestream.")
	    .add_property("n_samples", &G3Timestream_nsamples,
	      "Number of samples in the timestream. Equivalent to len(ts)")
	    .def("_assert_congruence", G3Timestream_assert_congruence,
	      "log_fatal() if units, length, start, or stop times do not match")
	    .def("_cxxslice", G3Timestream_getslice, "Slice-only __getitem__")
	    // Operators bound in python through numpy
	;
	scitbx::boost_python::container_conversions::from_python_sequence<G3Timestream, scitbx::boost_python::container_conversions::variable_capacity_policy>();
	register_pointer_conversions<G3Timestream>();

	// Add buffer protocol interface
	PyTypeObject *tsclass = (PyTypeObject *)ts.ptr();
	timestream_bufferprocs.bf_getbuffer = G3Timestream_getbuffer;
	tsclass->tp_as_buffer = &timestream_bufferprocs;
#if PY_MAJOR_VERSION < 3
	tsclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	bp::object tsm =
	  EXPORT_FRAMEOBJECT(G3TimestreamMap, init<>(), "Collection of timestreams indexed by logical detector ID")
	    .def(bp::std_map_indexing_suite<G3TimestreamMap, true>())
	    .def("CheckAlignment", &G3TimestreamMap::CheckAlignment)
	    .add_property("start", &G3TimestreamMap::GetStartTime,
	      &timestream_map_set_start_time,
	      "Time of the first sample in the time stream")
	    .add_property("stop", &G3TimestreamMap::GetStopTime,
	      &timestream_map_set_stop_time,
	      "Time of the final sample in the time stream")
	    .add_property("sample_rate", &G3TimestreamMap::GetSampleRate,
	      "Computed sample rate of the timestream.")
	    .add_property("n_samples", &G3TimestreamMap::NSamples,
	      "Number of samples in the timestream. Equivalent to the length "
	      "of one of the timestreams.")
	;
	register_pointer_conversions<G3TimestreamMap>();

	// Add buffer protocol interface
	PyTypeObject *tsmclass = (PyTypeObject *)tsm.ptr();
	timestreammap_bufferprocs.bf_getbuffer = G3TimestreamMap_getbuffer;
	timestreammap_bufferprocs.bf_releasebuffer = G3TimestreamMap_relbuffer;
	tsmclass->tp_as_buffer = &timestreammap_bufferprocs;
#if PY_MAJOR_VERSION < 3
	tsmclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif
}

