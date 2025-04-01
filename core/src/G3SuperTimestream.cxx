#include <serialization.h>
#include "G3SuperTimestream.h"

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

template <class A> void G3SuperTimestream::load(A &ar, unsigned v)
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

	G3Time start = times[0];
	G3Time stop = times[times.size() - 1];
	// XXX check times for missing samples

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
	size_t size = ndet * nsamp;

	bool is_floaty = (type_num == TYPE_NUM_FLOAT32 || type_num == TYPE_NUM_FLOAT64);

	ar & make_nvp("data_algo", data_algo);
	char *buf;

	// XXX handle nans and missing samples
	switch (type_num) {
	case TYPE_NUM_INT32: {
		std::shared_ptr<int32_t[]> data(new int32_t[size]);
		FromBuffer(names, nsamp, data, start, stop);
		buf = (char *)data.get();
		break;
	}
	case TYPE_NUM_FLOAT32: {
		std::shared_ptr<float[]> data(new float[size]);
		FromBuffer(names, nsamp, data, start, stop);
		buf = (char *)data.get();
		break;
	}
	case TYPE_NUM_INT64: {
		std::shared_ptr<int64_t[]> data(new int64_t[size]);
		FromBuffer(names, nsamp, data, start, stop);
		buf = (char *)data.get();
		break;
	}
	case TYPE_NUM_FLOAT64: {
		std::shared_ptr<double[]> data(new double[size]);
		FromBuffer(names, nsamp, data, start, stop);
		buf = (char *)data.get();
		break;
	}
	default:
		log_fatal("Invalid type num %d", (int) type_num);
	}

	SetFLACCompression((type_num == TYPE_NUM_INT32) ? flac_level : 0);

	if (data_algo == ALGO_NONE) {
		ar & make_nvp("data_raw", binary_data(buf, nbytes));
		return;
	}

	int count;
	std::vector<int> offsets;
	std::vector<double> quanta;

	// Read the flacblock
	ar & make_nvp("quanta", quanta);
	ar & make_nvp("offsets", offsets);
	ar & make_nvp("payload_bytes", count);
	std::unique_ptr<char []> cbuf(new char[count]);
	ar & make_nvp("payload", binary_data(&cbuf[0], count));

	// Decompress or copy into a buffer.
	FLAC__StreamDecoderWriteCallback this_write_callback;
	void (*expand_func)(struct flac_helper *, char*) = nullptr;
	void (*broadcast_func)(struct flac_helper *) = nullptr;
	void (*rescale_func)(struct flac_helper *, double) = nullptr;

	if (!is_64bit) {
		this_write_callback = &write_callback_int<int32_t>;
		expand_func = expand_branch<int32_t>;
		broadcast_func = broadcast_val<int32_t>;
		if (is_floaty)
			rescale_func = rescale<int32_t, float>;
	} else {
		this_write_callback = &write_callback_int<int64_t>;
		expand_func = expand_branch<int64_t>;
		broadcast_func = broadcast_val<int64_t>;
		if (is_floaty)
			rescale_func = rescale<int64_t, double>;
	}

#pragma omp parallel
	{
		// Each OMP thread needs its own workspace, FLAC decoder, and helper structure
		std::unique_ptr<char []> temp(new char[nsamp * elsize + 1]);
		FLAC__StreamDecoder *decoder = nullptr;
		struct flac_helper helper;

#pragma omp for
		for (int i=0; i<ndet; i++) {
			char* this_data = buf + elsize * nsamp * i;

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
		}

		if (decoder != nullptr)
			FLAC__stream_decoder_delete(decoder);

	} // omp parallel
}

#else

template <class A> void G3SuperTimestream::load(A &ar, unsigned v)
{
	log_fatal("Library missing FLAC or BZip2 compression");
}

#endif

// Save directly to a G3TimestreamMap object
template <class A> void G3SuperTimestream::save(A &ar, unsigned v) const
{
	log_fatal("Convert to G3Timestream map to serialize");
}

G3_SPLIT_SERIALIZABLE_CODE(G3SuperTimestream);
