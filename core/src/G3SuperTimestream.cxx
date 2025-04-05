#include <G3Logging.h>
#include <G3Timestream.h>
#include <serialization.h>

/*
  This is a minimal implementation of the Simons Observatory timestream storage
  class, to enable reading legacy data from disk and convert directly to
  G3TimestreamMap.  This class is only visible to the serialization library.
*/

class G3SuperTimestream : public G3TimestreamMap {
private:
	G3SuperTimestream() {};

	friend class cereal::access;

	template <class A> void load(A &ar, unsigned v) {
		log_fatal("Library missing FLAC or BZip2 compression");
	}

	template <class A> void save(A &ar, unsigned v) const {
		log_fatal("Convert to G3TimestreamMap to serialize");
	}

	SET_LOGGER("G3SuperTimestream");
};

namespace cereal {
	template <class A> struct specialize<A, G3SuperTimestream,
	    cereal::specialization::member_load_save> {};

	// Convert to G3TimestreamMap before handing off to the user
	namespace detail {
		template <> inline std::shared_ptr<void>
		PolymorphicCasters::upcast(std::shared_ptr<G3SuperTimestream> const & dptr,
		    std::type_info const & baseInfo) {
			return PolymorphicCasters::upcast(
			    std::make_shared<G3TimestreamMap>(*dptr), baseInfo);
		}
	}
}

G3_SERIALIZABLE(G3SuperTimestream, 0);

// Everything below is just the deserialization implementation, largely copied
// wholesale from the simonsobs/so3g library, and rearranged to handle decompression
// immediately on load.

#if defined(G3_HAS_FLAC) && defined(BZIP2_FOUND)

#ifdef _OPENMP
#include <omp.h>
#endif
#include <FLAC/stream_encoder.h>
#include <bzlib.h>

// Hard-code numpy type enums
#define TYPE_NUM_INT32 5
#ifdef __LP64__
#define TYPE_NUM_INT64 7
#else
#define TYPE_NUM_INT64 9
#endif
#define TYPE_NUM_FLOAT32 11
#define TYPE_NUM_FLOAT64 12

enum algos {
	ALGO_NONE = 0,
	ALGO_DO_FLAC = (1 << 0),
	ALGO_DO_BZ = (1 << 1),
	ALGO_DO_CONST = (1 << 2)
};

struct flac_helper {
	int bytes_remaining;
	char *src;
	char *dest;
	int start;
	int count;
};

FLAC__StreamDecoderReadStatus read_callback(const FLAC__StreamDecoder *decoder,
    FLAC__byte buffer[], size_t *bytes, void *client_data)
{
	auto fh = (struct flac_helper *)client_data;
	/* printf(" ... read %i (remaining: %i)\n", *bytes, fh->bytes_remaining); */
	if (fh->bytes_remaining == 0) {
		*bytes = 0;
		return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM;
	}
	if ((size_t)fh->bytes_remaining < *bytes)
		*bytes = fh->bytes_remaining;
	memcpy(buffer, fh->src, *bytes);
	fh->bytes_remaining -= *bytes;
	fh->src += *bytes;
	return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
}

template <typename T>
FLAC__StreamDecoderWriteStatus write_callback_int(const FLAC__StreamDecoder *decoder,
    const FLAC__Frame *frame, const FLAC__int32 *const buffer[], void *client_data)
{
	auto fh = (struct flac_helper *)client_data;
	int n = frame->header.blocksize;
	int drop = fh->start;
	if (drop >= n) {
		fh->start -= n;
	} else {
		n -= drop;
		fh->start = 0;
		if (n > fh->count)
			n = fh->count;
		for (int i=0; i<n; i++)
			((T*)fh->dest)[i] = buffer[0][i+drop];
		fh->dest += n * sizeof(T);
		fh->count -= n;
	}
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

static
void bz2_error_cb(int err) {
	switch(err) {
	case BZ_CONFIG_ERROR:
		log_fatal("BZ_CONFIG_ERROR (library compilation issue)");
	case BZ_PARAM_ERROR:
		log_fatal("BZ_PARAM_ERROR (bad blocksize, verbosity, etc)");
	case BZ_MEM_ERROR:
		log_fatal("BZ_MEM_ERROR (not enough memory is available)");
	case BZ_OUTBUFF_FULL:
		log_fatal("BZ_OUTBUFF_FULL (compressed data too long for buffer)");
	case BZ_OK:
		return;
	default:
		log_fatal("Unknown BZ error code %d", err);
	}
}

template <typename T>
void expand_branch(struct flac_helper *fh, char *temp)
{
	int n_bytes = fh->count * sizeof(T);
	unsigned int temp_size = n_bytes;

	int err = BZ2_bzBuffToBuffDecompress(temp, &temp_size, fh->src, n_bytes, 1, 0);
	if (err != BZ_OK)
		bz2_error_cb(err);
	// Add it in ...
	for (int i=0; i<fh->count; i++)
		((T*)fh->dest)[i] += ((T*)temp)[i + fh->start];
}

template <typename T>
void broadcast_val(struct flac_helper *fh)
{
	T val = *((T*)(fh->src));
	T *dest = (T*)fh->dest;
	for (int i=0; i<fh->count; i++)
		dest[i] += val;
}

inline int32_t _read_size(struct flac_helper *fh)
{
	auto v = *((int32_t*)(fh->src));
	fh->src += sizeof(v);
	return v;
}

template <typename Tin, typename Tout>
void rescale(struct flac_helper *fh, double scale)
{
	auto src = (Tin*)fh->dest;
	auto dest = (Tout*)fh->dest;
	for (int j=0; j<fh->count; j++)
		dest[j] = (Tout)src[j] * scale;
}

std::vector<bool> find_time_gaps(const G3VectorTime &times)
{
	std::vector<double> dts;
	dts.reserve(times.size() - 1);
	for (size_t i = 1; i < times.size(); i++)
		dts.push_back(times[i].time - times[i - 1].time);

	double dt;
	{
		int mid = dts.size() / 2;
		std::vector<double> sorted = dts;
		std::nth_element(sorted.begin(), sorted.begin() + mid, sorted.end());
		dt = sorted[mid];
	}

	std::vector<bool> gaps;
	gaps.reserve(times.size());
	gaps.push_back(false);

	for (size_t i = 0; i < dts.size(); i++) {
		int missing = (int)(std::round(dts[i] / dt - 1));
		for (int j = 0; j < missing; j++)
			gaps.push_back(true);
		gaps.push_back(false);
	}

	return gaps;
}

template <typename T>
void fill_gaps(struct flac_helper *fh, const std::vector<bool> &gaps, double fillval)
{
	std::vector<T> src((T*)fh->dest, (T*)fh->dest + fh->count);
	T *dest = (T*)fh->dest;

	int sidx = 0;
	for (size_t didx = 0; didx < gaps.size(); didx++)
		dest[didx] = (!gaps[didx] && sidx < fh->count) ? src[sidx++] : fillval;
}

template <>
void G3SuperTimestream::load(cereal::PortableBinaryInputArchive &ar, unsigned v)
{
	G3_CHECK_VERSION(v);
	using namespace cereal;

	ar & make_nvp("parent", base_class<G3FrameObject>(this));

	int8_t flac_level, bz2_workFactor, times_algo, data_algo;
	G3VectorTime times;
	G3VectorString names;

	// Compression options.
	ar & make_nvp("flac_level", flac_level);
	ar & make_nvp("bz2_workFactor", bz2_workFactor);

	ar & make_nvp("times_algo", times_algo);
	if (times_algo == ALGO_NONE) {
		ar & make_nvp("times", times);
	} else if (times_algo == ALGO_DO_BZ) {
		int n_samps;
		unsigned int max_bytes;
		ar & make_nvp("n_samps", n_samps);
		ar & make_nvp("comp_bytes", max_bytes);

		std::unique_ptr<char []> buf(new char[max_bytes]);
		ar & make_nvp("times_data", binary_data(&buf[0], max_bytes));

		std::vector<int64_t> ints(n_samps);
		unsigned int n_decomp = n_samps * sizeof(ints[0]);
		int err = BZ2_bzBuffToBuffDecompress((char*)&ints[0], &n_decomp, &buf[0],
		    max_bytes, 1, 0);
		if (err != BZ_OK)
			bz2_error_cb(err);
		times = G3VectorTime(ints.begin(), ints.end());
	} else
		log_fatal("No support for compression algorithm times_algo=%d", (int)times_algo);

	// Create a mask of missing samples
	G3Time start = times[0];
	G3Time stop = times[times.size() - 1];
	std::vector<bool> time_gaps = find_time_gaps(times);

	ar & make_nvp("names", names);

	// Read the desc.
	intptr_t type_num, ndim, shape[32], nbytes;
	ar & make_nvp("type_num", type_num);
	ar & make_nvp("ndim", ndim);
	ar & make_nvp("shape", shape);
	ar & make_nvp("nbytes", nbytes);

	bool is_64bit = (type_num == TYPE_NUM_INT64 || type_num == TYPE_NUM_FLOAT64);
	size_t elsize = is_64bit ? sizeof(int64_t) : sizeof(int32_t);
	int ndet = shape[0];
	int nsamp = shape[1];
	int nsamp_filled = time_gaps.size();
	size_t size = ndet * nsamp_filled;

	// Add bad sample mask to map as an extra detector channel at the end
	bool is_floaty = (type_num == TYPE_NUM_FLOAT32 || type_num == TYPE_NUM_FLOAT64);
	bool found_gaps = (nsamp_filled != nsamp);
	size_t t0 = size;
	if (found_gaps && !is_floaty) {
		names.push_back("_nanmask");
		size += nsamp_filled;
	}

	ar & make_nvp("data_algo", data_algo);
	char *buf;

	switch (type_num) {
	case TYPE_NUM_INT32: {
		std::shared_ptr<int32_t[]> data(new int32_t[size]);
		FromBuffer(names, nsamp_filled, data, start, stop);
		buf = (char *)data.get();

		if (found_gaps)
			for (size_t j = 0, i = t0; i < size; i++, j++)
				data[i] = (int32_t)time_gaps[j];
		break;
	}
	case TYPE_NUM_FLOAT32: {
		std::shared_ptr<float[]> data(new float[size]);
		FromBuffer(names, nsamp_filled, data, start, stop);
		buf = (char *)data.get();
		break;
	}
	case TYPE_NUM_INT64: {
		std::shared_ptr<int64_t[]> data(new int64_t[size]);
		FromBuffer(names, nsamp_filled, data, start, stop);
		buf = (char *)data.get();

		if (found_gaps)
			for (size_t j = 0, i = t0; i < size; i++, j++)
				data[i] = (int64_t)time_gaps[j];
		break;
	}
	case TYPE_NUM_FLOAT64: {
		std::shared_ptr<double[]> data(new double[size]);
		FromBuffer(names, nsamp_filled, data, start, stop);
		buf = (char *)data.get();
		break;
	}
	default:
		log_fatal("Invalid type num %d", (int) type_num);
	}

	SetFLACCompression((type_num == TYPE_NUM_INT32) ? flac_level : 0);

	if (data_algo == ALGO_NONE && !found_gaps) {
		ar & make_nvp("data_raw", binary_data(buf, nbytes));
		return;
	}

	int count = nbytes;
	std::vector<int> offsets;
	std::vector<double> quanta;

	if (data_algo == ALGO_NONE) {
		// fill gaps in timestreams in parallel
		for (int i=0; i<ndet; i++)
			offsets.push_back(elsize * nsamp * i);
	} else {
		// Read the flac block
		ar & make_nvp("quanta", quanta);
		ar & make_nvp("offsets", offsets);
		ar & make_nvp("payload_bytes", count);
	}

	std::unique_ptr<char []> cbuf(new char[count]);
	ar & make_nvp("payload", binary_data(&cbuf[0], count));

	// Decompress or copy into a buffer.
	FLAC__StreamDecoderWriteCallback this_write_callback;
	void (*expand_func)(struct flac_helper *, char*) = nullptr;
	void (*broadcast_func)(struct flac_helper *) = nullptr;
	void (*rescale_func)(struct flac_helper *, double) = nullptr;
	void (*gapfill_func)(struct flac_helper *, const std::vector<bool> &, double) = nullptr;

	if (!is_64bit) {
		this_write_callback = &write_callback_int<int32_t>;
		expand_func = expand_branch<int32_t>;
		broadcast_func = broadcast_val<int32_t>;
		if (is_floaty)
			rescale_func = rescale<int32_t, float>;
		gapfill_func = is_floaty ? fill_gaps<float> : fill_gaps<int32_t>;
	} else {
		this_write_callback = &write_callback_int<int64_t>;
		expand_func = expand_branch<int64_t>;
		broadcast_func = broadcast_val<int64_t>;
		if (is_floaty)
			rescale_func = rescale<int64_t, double>;
		gapfill_func = is_floaty ? fill_gaps<double> : fill_gaps<int64_t>;
	}

#pragma omp parallel
	{
		// Each OMP thread needs its own workspace, FLAC decoder, and helper structure
		std::unique_ptr<char []> temp(new char[nsamp * elsize + 1]);
		FLAC__StreamDecoder *decoder = nullptr;
		struct flac_helper helper;

#pragma omp for
		for (int i=0; i<ndet; i++) {
			char* this_data = buf + elsize * nsamp_filled * i;

			// Cue up this detector's data and read the algo code.
			helper.src = &cbuf[0] + offsets[i];
			int8_t algo = (data_algo == ALGO_NONE) ? ALGO_NONE : *(helper.src++);

			helper.dest = this_data;
			helper.start = 0;
			helper.count = nsamp;

			if (algo == ALGO_NONE)
				memcpy(helper.dest, helper.src, nsamp * elsize);

			if (algo & ALGO_DO_FLAC) {
				if (decoder == nullptr)
					decoder = FLAC__stream_decoder_new();
				helper.bytes_remaining = _read_size(&helper);

				FLAC__stream_decoder_init_stream(decoder, read_callback,
				    NULL, NULL, NULL, NULL, *this_write_callback, NULL,
				    flac_decoder_error_cb, (void*)&helper);
				FLAC__stream_decoder_process_until_end_of_stream(decoder);
				FLAC__stream_decoder_finish(decoder);

				// reset output pointers
				helper.dest = this_data;
				helper.start = 0;
				helper.count = nsamp;
			}

			// A bz2 field of slow offsets?
			if (algo & ALGO_DO_BZ) {
				helper.bytes_remaining = _read_size(&helper);
				expand_func(&helper, &temp[0]);
			}

			// Single flat offset?
			if (algo & ALGO_DO_CONST)
				broadcast_func(&helper);

			// Now convert for precision.
			if (data_algo != ALGO_NONE && is_floaty)
				rescale_func(&helper, quanta[i]);

			// Fill missing samples with zeros or nans
			if (found_gaps)
				gapfill_func(&helper, time_gaps, is_floaty ? NAN : 0.0);
		}

		if (decoder != nullptr)
			FLAC__stream_decoder_delete(decoder);

	} // omp parallel
}

#endif

