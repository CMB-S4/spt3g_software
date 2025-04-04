//#define NO_IMPORT_ARRAY

#include <pybindings.h>
#include <iostream>
#include <boost/python.hpp>

#include <container_pybindings.h>
#include <so/G3SuperTimestream.h>
#include <cereal/types/utility.hpp>

#ifdef _OPENMP
# include <omp.h>
#endif // ifdef _OPENMP
#include <FLAC/stream_encoder.h>
#include <bzlib.h>

// Debugging variables for compressors.
// BZ2_VERBOSITY can range from 0 (silent) to 4.
#define SO3G_BZ2_VERBOSITY 0
#define SO3G_BZ2_BLOCKSIZE 5

static
std::string get_bz2_error_string(int err) {
    std::ostringstream s;
    switch(err) {
    case BZ_CONFIG_ERROR:
        s << "BZ_CONFIG_ERROR (library compilation issue)";
        break;
    case BZ_PARAM_ERROR:
        s <<  "BZ_PARAM_ERROR (bad blocksize, verbosity, etc)";
        break;
    case BZ_MEM_ERROR:
        s << "BZ_MEM_ERROR (not enough memory is available)";
        break;
    case BZ_OUTBUFF_FULL:
        s << "BZ_OUTBUFF_FULL (compressed data too long for buffer)";
        break;
    case BZ_OK:
        s << "BZ_OK (no problem)";
        break;
    default:
        s << "Unknown BZ error code " << err;
    }
    return s.str();
}

static
std::string get_algo_error_string(std::string var_name, int algo_code)
{
    std::ostringstream s;
    s << "No support for compression algorithm " << var_name << "=" << algo_code;
    return s.str();
}

// Split each datum in x into:
//
//    x[i] = r[i] + y[i]
//
// where r is a multiple of snap_to, -step_at <= y[i] < step_at.  The
// algorithm tries to have r[i] change slowly, meaning that r[i+1]
// will not differ from r[i] unless necessary to satisfy the
// range restriction on y[i].
//
// This is only implemented for snap_to and step_at both powers of 2!
//
// Let's allow r and x to point to same memory.
//
// This will return the number of values that y takes; 0 if there are
// no data, 1 if it's a constant value (even if that value is 0), or
// more than 1 for more than 1.
template <typename T>
static
int rebranch(int32_t *y, T *r, T *x, int n_samps, T snap_to, T step_at)
{
    T branch = 0;
    int branch_count = 0;
    for (int i=0; i < n_samps; i++) {
        T new_branch = (x[i] + snap_to / 2) & ~(snap_to - 1);
        if (
            new_branch - branch >= step_at ||
            new_branch - branch <  -step_at ||
            branch_count == 0
        ) {
            branch = new_branch;
            branch_count++;
        }
        y[i] = x[i] - branch;
        // Don't update r until you use x -- they are allowed to overlap.
        r[i] = branch;
    }
    return branch_count;
}

FLAC__StreamEncoderWriteStatus flac_encoder_write_cb(
    const FLAC__StreamEncoder *encoder,
    const FLAC__byte buffer[],
    size_t bytes,
    unsigned samples,
    unsigned current_frame,
    void *client_data)
{
    auto ablob = (struct G3SuperTimestream::array_blob *)client_data;
    if ((size_t)ablob->count + bytes > (size_t)ablob->size)
        return FLAC__STREAM_ENCODER_WRITE_STATUS_FATAL_ERROR;
    memcpy(ablob->buf + ablob->count, buffer, bytes);
    ablob->count += bytes;
    return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

inline int32_t* reserve_size_field(struct G3SuperTimestream::array_blob *ablob)
{
    if ((size_t)(ablob->size - ablob->count) < sizeof(int32_t))
        return nullptr;
    auto p = (int32_t*)(ablob->buf + ablob->count);
    ablob->count += sizeof(int32_t);
    *p = ablob->count;
    return p;
}
inline bool close_size_field(struct G3SuperTimestream::array_blob *ablob, int32_t *dest)
{
    if (dest == NULL)
        return false;
    *dest = ablob->count - *dest;
    return true;
}

static
struct G3SuperTimestream::array_blob encode_array(
    PyArrayObject *array, std::vector<double> quanta,
    G3SuperTimestream::options_type options)
{
    struct G3SuperTimestream::array_blob ablob;

    // We will write (possibly) compressed data for n_det
    // detectors by n samples into a buffer.  The data block for
    // detector i starts at offsets[i] and ends at offsets[i+1].
    //
    // Each data block has the structure:
    //   [data_block] = [algo_code] [encoded_block1] [encoded_block2]
    //
    // The algo_code is a single byte drawing values from the
    // algos enum/bitmask that describes what encoded blocks will
    // follow.
    //
    // If algo_code=0 (a.k.a. ALGO_NONE) then there is one
    // encoded_block and it is simply the flat uncompressed binary
    // data for that channel.
    //
    // If algo_code != 0, then the encoded blocks will consist
    // first of a FLAC block (if algo_code & ALGO_DO_FLAC),
    // followed by a BZ2 block (if ALGO_DO_BZ) or a CONST block
    // (if ALGO_DO_CONST).
    //
    // The FLAC block has the format [length] [flac_data].  The
    // length is an int32_t giving the number of bytes in
    // flac_data.  The flac_data is the encoding of n samples of
    // single-channel 24-bit data.
    //
    // The BZ2 block has the format [length] [bz2_data].  The
    // length is an int32_t giving the number of bytes in
    // bz2_data.  The bz2_data is the encoding of n samples of the
    // full-width data (i.e. n*sizeof(dtype) bytes).
    //
    // The CONST block is a single full-width value, [datum].
    // This value represents an offset to add to the data.
    //
    // The FLAC data, if present, decode to int32_t.  The BZ2
    // data, if present, decode to int32_t or int64_t.  The CONST
    // datum, if present, represents a single int32_t or int64_t.
    // These two vectors and single offset are added together to
    // form the output integer array.
    //
    // (In practice we will only have one or the other of BZ2 and
    // CONST blocks.  In the common case that CONST would decode
    // to 0, it will not be included at all.)

    int n_chans = PyArray_SHAPE(array)[0];
    int n_samps = PyArray_SHAPE(array)[1];
    int itemsize = PyArray_ITEMSIZE(array);

    // Max bytes needed to store all the "compressed" data.
    int n_max = PyArray_NBYTES(array) + n_chans;

    // Initialize the buffer.
    ablob.buf = new char[n_max];
    ablob.size = n_max;
    ablob.count = 0;
    ablob.offsets.push_back(0);

    int32_t d[n_samps];
    const int32_t *chan_ptrs[1] = {d};

    char r[n_samps * sizeof(int64_t)];
    auto r32 = (int32_t*)r;
    auto r64 = (int64_t*)r;

    npy_intp type_num = PyArray_TYPE(array);
    char *src = (char*)(PyArray_DATA(array));

    int32_t M = (1 << 24);
    for (int i=0; i<n_chans; i++, src+=PyArray_STRIDES(array)[0]) {
        // Save the current ablob.buf pointer in case we need to
        // roll back the compression, which we will do if we
        // end up writing to block_limit or beyond.
        int block_start = ablob.count++;
        int block_limit = ablob.count + n_samps * itemsize;

        // Oh, the things we'll try.
        bool try_flac = (options.data_algo & G3SuperTimestream::ALGO_DO_FLAC);
        bool try_bz = (options.data_algo & G3SuperTimestream::ALGO_DO_BZ);
        bool try_const = false;

        int8_t algo_used = 0;
        bool ok = true;

        if (try_flac) {
            int32_t *flac_size = reserve_size_field(&ablob);
            int branches = 0;
            if (type_num == NPY_INT32) {
                branches = rebranch<int32_t>(d, r32, (int32_t*)src, n_samps, M, M/2);
            } else if (type_num == NPY_INT64) {
                branches = rebranch<int64_t>(d, r64, (int64_t*)src, n_samps, M, M/2);
            } else if (type_num == NPY_FLOAT32) {
                for (int j=0; j < n_samps; j++)
                    r32[j] = roundf(((float*)src)[j] / quanta[i]);
                branches = rebranch<int32_t>(d, r32, r32, n_samps, M, M/2);
            } else if (type_num == NPY_FLOAT64) {
                for (int j=0; j < n_samps; j++)
                    r64[j] = round(((double*)src)[j] / quanta[i]);
                branches = rebranch<int64_t>(d, r64, r64, n_samps, M, M/2);
            } else
                throw g3supertimestream_exception("Invalid array type encountered.");

            // Re-using this encoder for all
            // detectors... seems to not work if process
            // or finish exit with failure.
            auto encoder = FLAC__stream_encoder_new();
            FLAC__stream_encoder_set_channels(encoder, 1);
            FLAC__stream_encoder_set_bits_per_sample(encoder, 24);
            FLAC__stream_encoder_set_compression_level(encoder, options.flac_level);
            FLAC__stream_encoder_init_stream(
                encoder, flac_encoder_write_cb, NULL, NULL, NULL, (void*)(&ablob));

            // If encoding fails, check that it looks like
            // a simple buffer-to-small error, then continue.
            if (!FLAC__stream_encoder_process(encoder, chan_ptrs, n_samps) ||
                !FLAC__stream_encoder_finish(encoder)) {
                auto state = FLAC__stream_encoder_get_state(encoder);
                if (state != FLAC__STREAM_ENCODER_CLIENT_ERROR) {
                    const char *estr = FLAC__stream_encoder_get_resolved_state_string(encoder);
                    std::ostringstream s;
                    s << "Unexpected FLAC encoder error: " << estr
                      << " on channel " << i << " (" << n_samps << " samples)";
                    throw g3supertimestream_exception(s.str());
                }
                // Discard any partial encoding ... let bzip have a go.
                ablob.count = block_start + 1;
                        } else {
                // Great, update FLAC block size and record it.
                ok = close_size_field(&ablob, flac_size);
                algo_used |= G3SuperTimestream::ALGO_DO_FLAC;

                // Const remainder?
                if (branches <= 1) {
                    try_bz = false;
                    try_const = true;
                }
            }
            FLAC__stream_encoder_delete(encoder);
        }

        if (!(algo_used & G3SuperTimestream::ALGO_DO_FLAC)) {
            // Just copy the raw data into the r buffer,
            // where the bz code will pick it up.
            memcpy(r, src, n_samps * itemsize);
            if (type_num == NPY_FLOAT32) {
                for (int j=0; j<n_samps; j++)
                    ((int32_t*)r)[j] = roundf(((float*)r)[j] / quanta[i]);
            } else if (type_num == NPY_FLOAT64) {
                for (int j=0; j<n_samps; j++)
                    ((int64_t*)r)[j] = round(((double*)r)[j] / quanta[i]);
            }
        }

        if (ok && try_const) {
            // Note "try_const" means "the data are
            // constant, it's ok to treat them as such".
            // The code below doesn't check the
            // assumption.
            //
            // Take care in the cute handling of special
            // case r[0] == 0.  If you naively leave out
            // the ALGO_DO_CONST bit, and if FLAC
            // compression was disabled, then you end up
            // storing algo=0, which the decoder will
            // interpret as raw format, not zero-byte
            // ultracompression representing all zeros!
            int n_copy = 0;
            switch (type_num) {
            case NPY_FLOAT32:
            case NPY_INT32:
                if (r32[0] != 0 || algo_used == 0)
                    n_copy = sizeof(int32_t);
                break;
            case NPY_FLOAT64:
            case NPY_INT64:
                if (r64[0] != 0 || algo_used == 0)
                    n_copy = sizeof(int64_t);
                break;
            }
            if (ablob.count + n_copy < block_limit) {
                if (n_copy > 0) {
                    memcpy(ablob.buf + ablob.count, r, n_copy);
                    ablob.count += n_copy;
                    algo_used |= G3SuperTimestream::ALGO_DO_CONST;
                }
            } else {
                ok = false;
            }
        }
        if (ok && try_bz) {
            int32_t *bz_size = reserve_size_field(&ablob);
            algo_used |= G3SuperTimestream::ALGO_DO_BZ;
            if (ablob.count < block_limit) {
                unsigned int n_write = block_limit - ablob.count;
                int err = BZ2_bzBuffToBuffCompress(
                    ablob.buf + ablob.count, &n_write, r, n_samps * itemsize,
                    SO3G_BZ2_BLOCKSIZE, SO3G_BZ2_VERBOSITY, options.bz2_workFactor);
                if (err == BZ_OUTBUFF_FULL) {
                    // Too long, don't bother.
                    ok = false;
                } else if (err != BZ_OK) {
                    throw g3supertimestream_exception(get_bz2_error_string(err));
                } else {
                    ablob.count += n_write;
                    ok = close_size_field(&ablob, bz_size);
                }
            } else
                ok = false;
        }

        if (ok && ablob.count < block_limit) {
            // That went well.
            ablob.buf[block_start] = algo_used;
        } else {
            // Bail out into raw copy.
            ablob.buf[block_start] = G3SuperTimestream::ALGO_NONE;
            memcpy(ablob.buf + block_start + 1, src, n_samps * itemsize);
            ablob.count = block_limit;
        }
        ablob.offsets.push_back(ablob.count);
    }

    return ablob;
}


/* G3SuperTimestream */

std::string G3SuperTimestream::Description() const
{
    std::ostringstream s;
    s << "G3SuperTimestream("
      << names.size() << ", " << times.size() << ")";
    return s.str();
}

std::string G3SuperTimestream::Summary() const
{
    return Description();
}

template <class A> void G3SuperTimestream::load(A &ar, unsigned v)
{
    G3_CHECK_VERSION(v);
    using namespace cereal;
    ar & make_nvp("parent", base_class<G3FrameObject>(this));

        // Compression options.
    ar & make_nvp("flac_level", options.flac_level);
    ar & make_nvp("bz2_workFactor", options.bz2_workFactor);

    ar & make_nvp("times_algo", options.times_algo);
    if (options.times_algo == ALGO_NONE) {
        ar & make_nvp("times", times);
    } else if (options.times_algo == ALGO_DO_BZ) {
        int n_samps;
        unsigned int max_bytes;
        ar & make_nvp("n_samps", n_samps);
        ar & make_nvp("comp_bytes", max_bytes);

                char _buf[max_bytes];
                char *buf = (char*)_buf;
        ar & make_nvp("times_data", binary_data(buf, max_bytes));

        std::vector<int64_t> ints(n_samps);
        unsigned int n_decomp = n_samps * sizeof(ints[0]);
        int err = BZ2_bzBuffToBuffDecompress(
                    (char*)&ints[0], &n_decomp, buf, max_bytes,
                 1, SO3G_BZ2_VERBOSITY);
        if (err != BZ_OK)
            throw g3supertimestream_exception(get_bz2_error_string(err));
        times = G3VectorTime(ints.begin(), ints.end());
    } else
        throw g3supertimestream_exception(
            get_algo_error_string("times_algo", options.times_algo));

    ar & make_nvp("names", names);

    // Read the desc.
        ar & make_nvp("type_num", desc.type_num);
        ar & make_nvp("ndim", desc.ndim);
    ar & make_nvp("shape", desc.shape);
        ar & make_nvp("nbytes", desc.nbytes);

    ar & make_nvp("data_algo", options.data_algo);
    if (options.data_algo == ALGO_NONE) {
        array = (PyArrayObject*)
            PyArray_SimpleNew(desc.ndim, desc.shape, desc.type_num);
        ar & make_nvp("data_raw", binary_data((char*)PyArray_DATA(array),
                              PyArray_NBYTES(array)));
    } else {
        // Read the flacblock
        ablob = new struct array_blob;
        ar & make_nvp("quanta", quanta);
        ar & make_nvp("offsets", ablob->offsets);
        ar & make_nvp("payload_bytes", ablob->count);
        ablob->buf = new char[ablob->count];
        ar & make_nvp("payload", binary_data(ablob->buf, ablob->count));
    }
    dataful = true;
    float_mode = (desc.type_num == NPY_FLOAT32 ||
              desc.type_num == NPY_FLOAT64);
}

template <class A> void G3SuperTimestream::save(A &ar, unsigned v) const
{
    using namespace cereal;
    ar & make_nvp("parent", base_class<G3FrameObject>(this));

        // Compression options.
    ar & make_nvp("flac_level", options.flac_level);
    ar & make_nvp("bz2_workFactor", options.bz2_workFactor);

    if (options.times_algo == ALGO_DO_BZ && times.size() > 0) {
        // Try to bz2 compress.  Convert to a vector of int64_t first.
        // Note this doesn't run unless n_samps > 0.
        auto time_ints = vector<long>(times.begin(), times.end());
        int n_samps = time_ints.size();
        unsigned int max_bytes = n_samps * sizeof(time_ints[0]);

                char _buf[max_bytes];
        char *buf = _buf;

        int err = BZ2_bzBuffToBuffCompress(
            buf, &max_bytes, (char*)&time_ints[0], max_bytes,
            SO3G_BZ2_BLOCKSIZE, SO3G_BZ2_VERBOSITY, options.bz2_workFactor);
        if (err == BZ_OUTBUFF_FULL) {
            // Fallback to no-compression.
            ar & make_nvp("times_algo", (int8_t)ALGO_NONE);
            ar & make_nvp("times", times);
        } else if (err != BZ_OK) {
            throw g3supertimestream_exception(get_bz2_error_string(err));
        } else {
            ar & make_nvp("times_algo", (int8_t)ALGO_DO_BZ);
            ar & make_nvp("n_samps", n_samps);
            ar & make_nvp("comp_bytes", max_bytes);
            ar & make_nvp("times_data", binary_data(buf, max_bytes));
        }
    } else {
        ar & make_nvp("times_algo", (int8_t)ALGO_NONE);
        ar & make_nvp("times", times);
    }

    ar & make_nvp("names", names);

    // Write the desc.
        ar & make_nvp("type_num", desc.type_num);
        ar & make_nvp("ndim", desc.ndim);
    ar & make_nvp("shape", desc.shape);
        ar & make_nvp("nbytes", desc.nbytes);

    ar & make_nvp("data_algo", options.data_algo);
    if (options.data_algo == ALGO_NONE) {
        if (array == nullptr)
            throw g3supertimestream_exception(
                "Unexpected state: array is NULL.");
        // Check the endianness
        //auto this_descr = PyArray_DESCR(array);
        if (!PyArray_EquivByteorders(PyArray_DESCR(array)->byteorder, NPY_LITTLE))
            throw g3supertimestream_exception(
                "The byte_order of the data array is not acceptable.");
        // Might as well use numpy to repack it properly...
        PyArrayObject *contig = PyArray_GETCONTIGUOUS(array);
        ar & make_nvp("data_raw", binary_data((char*)PyArray_DATA(contig),
                              PyArray_NBYTES(contig)));
        Py_DECREF((PyObject*)contig);
    } else {
        struct array_blob *_ablob = ablob;
        if (_ablob == nullptr) {
            // Encode to a copy.
            _ablob = new struct array_blob;
            *_ablob = encode_array(array, quanta, options);
        }

        // Write the ablobblock
        ar & make_nvp("quanta", quanta);
        ar & make_nvp("offsets", _ablob->offsets);
        ar & make_nvp("payload_bytes", _ablob->count);
        ar & make_nvp("payload", binary_data(_ablob->buf, _ablob->count));

        if (_ablob != ablob) {
            delete _ablob->buf;
            delete _ablob;
        }
    }
}

bool G3SuperTimestream::Encode() {
    if (array == nullptr)
        return false;

    // Compress the array data.
    if (options.data_algo == ALGO_NONE)
        return false;
    else {
        ablob = new struct array_blob;
        *ablob = encode_array(array, quanta, options);
        Py_XDECREF(array);
        array = nullptr;
    }
    return true;
}

struct flac_helper {
    int bytes_remaining;
    char *src;
    char *dest;
    int start;
    int count;
};

FLAC__StreamDecoderReadStatus read_callback(
    const FLAC__StreamDecoder *decoder, FLAC__byte buffer[], size_t *bytes, void *client_data)
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
FLAC__StreamDecoderWriteStatus write_callback_int(
    const FLAC__StreamDecoder *decoder, const FLAC__Frame *frame, const FLAC__int32 *const buffer[], void *client_data)
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
        printf("FLAC decoding error (lost sync)");
    case FLAC__STREAM_DECODER_ERROR_STATUS_BAD_HEADER:
        printf("FLAC decoding error (bad header)");
    case FLAC__STREAM_DECODER_ERROR_STATUS_FRAME_CRC_MISMATCH:
        printf("FLAC decoding error (CRC mismatch)");
    case FLAC__STREAM_DECODER_ERROR_STATUS_UNPARSEABLE_STREAM:
        printf("FLAC decoding error (unparseable stream)");
    default:
        printf("FLAC decoding error (%d)", status);
    }
}

template <typename T>
void expand_branch(struct flac_helper *fh, int n_bytes, char *temp)
{
    unsigned int temp_size = n_bytes;

    int err = BZ2_bzBuffToBuffDecompress(
        temp, &temp_size, fh->src, n_bytes,
        1, SO3G_BZ2_VERBOSITY);
    if (err != BZ_OK)
        throw g3supertimestream_exception(get_bz2_error_string(err));
    // Add it in ...
    for (int i=0; i<fh->count; i++)
        ((T*)fh->dest)[i] += ((T*)temp)[i + fh->start];
}

template <typename T>
void broadcast_val(struct flac_helper *fh, int nsamps)
{
    T val = *((T*)(fh->src));
    T *dest = (T*)fh->dest;
    for (int i=0; i<nsamps; i++)
        dest[i] += val;
}

inline int32_t _read_size(struct flac_helper *fh)
{
    auto v = *((int32_t*)(fh->src));
    fh->src += sizeof(v);
    return v;
}

bool G3SuperTimestream::Decode()
{
    if (ablob == nullptr)
        return false;

    if (options.data_algo == ALGO_NONE)
        throw g3supertimestream_exception(
            "Decode called with ablob buffer but data_algo=0.");

    array = (PyArrayObject*)
        PyArray_ZEROS(desc.ndim, desc.shape, desc.type_num, 0);

    Extract(bp::object(bp::handle<>(bp::borrowed(reinterpret_cast<PyObject*>(array)))),
        bp::object(), bp::object(), 0, desc.shape[1]);

    // Destroy the flac bundle.
    delete ablob->buf;
    delete ablob;
    ablob = nullptr;

    return true;
}

bool G3SuperTimestream::Extract(
    bp::object dest, bp::object dest_indices, bp::object src_indices,
    int start, int stop)
{
    int n_det_ex = desc.shape[0];

    PyArrayObject *_dest = (PyArrayObject*)dest.ptr();

        auto _src_indices = BufferWrapper<int64_t>("src_indices", src_indices, true,
                                                   vector<int>{-1});
        if (_src_indices->obj != NULL)
            n_det_ex = _src_indices->shape[0];

        auto _dest_indices = BufferWrapper<int64_t>("dest_indices", dest_indices, true,
                                                   vector<int>{n_det_ex});

    // Sample index.
    if (start < 0 || start > desc.shape[1])
        throw g3supertimestream_exception(
                        "sample start index out of bounds");
    if (stop < 0)
        stop = desc.shape[1];
    if (stop < start || stop > desc.shape[1])
        throw g3supertimestream_exception(
                        "sample stop index out of bounds");

    int copy_count = stop - start;

    if (!PyArray_Check(_dest))
        throw g3supertimestream_exception(
            "Destination array must be ndarray.");
    if (PyArray_TYPE(_dest) != desc.type_num)
        throw g3supertimestream_exception(
            "Destination array not of correct type.");
    if ((PyArray_NDIM(_dest) != 2) ||
        (PyArray_SHAPE(_dest)[0] < n_det_ex) ||
        (PyArray_SHAPE(_dest)[1] != copy_count))
        throw g3supertimestream_exception(
            "Destination array does not have correct shape.");
    if (PyArray_STRIDE(_dest, 1) != PyArray_ITEMSIZE(_dest))
        throw g3supertimestream_exception(
            "Destination array should be strictly packed on last dimension.");

    // Insist array is not already decompressed.
    if (ablob == nullptr)
        throw g3supertimestream_exception(
            "These data are already fully decoded; use .data.");

    // Decompress or copy into a buffer.
    FLAC__StreamDecoderWriteCallback this_write_callback;
    void (*expand_func)(struct flac_helper *, int, char*) = nullptr;
    void (*broadcast_func)(struct flac_helper *, int) = nullptr;
    int elsize = 0;

    switch (desc.type_num) {
    case NPY_INT32:
    case NPY_FLOAT32:
        this_write_callback = &write_callback_int<int32_t>;
        expand_func = expand_branch<int32_t>;
        broadcast_func = broadcast_val<int32_t>;
        elsize = sizeof(int32_t);
        break;
    case NPY_INT64:
    case NPY_FLOAT64:
        this_write_callback = &write_callback_int<int64_t>;
        expand_func = expand_branch<int64_t>;
        broadcast_func = broadcast_val<int64_t>;
        elsize = sizeof(int64_t);
        break;
    default:
        throw g3supertimestream_exception("Invalid array type encountered.");
    }

#pragma omp parallel
    {

    // Each OMP thread needs its own workspace, FLAC decoder, and helper structure
    char temp[desc.shape[1] * elsize + 1];
    FLAC__StreamDecoder *decoder = nullptr;
    struct flac_helper helper;

#pragma omp for
    for (int i=0; i<n_det_ex; i++) {
        // Source vector index
        int src_i = i;
        if (_src_indices.test()) {
            src_i = *_src_indices.ptr_1d(i);
            if (src_i < 0 || src_i >= desc.shape[0])
                continue;
        }

        // Dest vector index
        int dest_i = i;
        if (_dest_indices.test()) {
            dest_i = *_dest_indices.ptr_1d(i);
            if (dest_i < 0 || dest_i >= PyArray_SHAPE(_dest)[0])
                continue;
        }

        char* this_data = (char*)PyArray_DATA(_dest) + PyArray_STRIDES(_dest)[0]*dest_i;

        // Cue up this detector's data and read the algo code.
        helper.src = ablob->buf + ablob->offsets[src_i];
        int8_t algo = *(helper.src++);

        if (algo == ALGO_NONE) {
            memcpy(this_data, helper.src + start * elsize, copy_count * elsize);
        }
        if (algo & ALGO_DO_FLAC) {
            if (decoder == nullptr)
                decoder = FLAC__stream_decoder_new();
            helper.bytes_remaining = _read_size(&helper);
            helper.dest = this_data;
            helper.start = start;
            helper.count = copy_count;

            FLAC__stream_decoder_init_stream(
                decoder, read_callback, NULL, NULL, NULL, NULL,
                *this_write_callback, NULL, flac_decoder_error_cb,
                (void*)&helper);
            FLAC__stream_decoder_process_until_end_of_stream(decoder);
            FLAC__stream_decoder_finish(decoder);
        }

        // A bz2 field of slow offsets?
        if (algo & ALGO_DO_BZ) {
            helper.bytes_remaining = _read_size(&helper);
            helper.dest = this_data;
            helper.start = start;
            helper.count = copy_count;
            expand_func(&helper, desc.shape[1] * elsize, (char*)temp);
        }

        // Single flat offset?
        if (algo & ALGO_DO_CONST) {
            helper.dest = this_data;
            broadcast_func(&helper, copy_count);
        }

        // Now convert for precision.
        if (desc.type_num == NPY_FLOAT32) {
            auto src = (int32_t*)this_data;
            auto dest = (float*)this_data;
            for (int j=0; j<PyArray_SHAPE(_dest)[1]; j++)
                dest[j] = (float)src[j] * quanta[src_i];
        } else if (desc.type_num == NPY_FLOAT64) {
            auto src = (int64_t*)this_data;
            auto dest = (double*)this_data;
            for (int j=0; j<PyArray_SHAPE(_dest)[1]; j++)
                dest[j] = src[j] * quanta[src_i];
        }
    }
    if (decoder != nullptr)
        FLAC__stream_decoder_delete(decoder);
    } // omp parallel

    return true;
}

int G3SuperTimestream::Options(
    int enable, int flac_level, int bz2_workFactor,
    int data_algo, int times_algo)
{
    // Global enable / disable first.
    if (enable == 0) {
        options.data_algo = 0;
        options.times_algo = 0;
    } else if (enable > 0) {
        options.data_algo = ALGO_DO_FLAC | ALGO_DO_BZ;
        options.times_algo = ALGO_DO_BZ;
    }

    if (data_algo >= 0)
        options.data_algo = data_algo;
    if (times_algo >= 0)
        options.times_algo = times_algo;

    if (flac_level >= 0)
        options.flac_level = flac_level;
    if (bz2_workFactor >= 0)
        options.bz2_workFactor = bz2_workFactor;

    return 0;
}


template <typename T>
static
void _apply_cals_typed(PyArrayObject *array, std::vector<double> cals)
{
    for (int i=0; i<PyArray_SHAPE(array)[0]; i++) {
        auto dest = (T*)((char*)PyArray_DATA(array) + PyArray_STRIDES(array)[0] * i);
        for (int j=0; j<PyArray_SHAPE(array)[1]; j++)
            dest[j] *= cals[i];
    }
}

void _apply_cals(PyArrayObject *array, std::vector<double> cals)
{
    switch (PyArray_TYPE(array)) {
    case NPY_FLOAT32:
        _apply_cals_typed<float>(array, cals);
        break;
    case NPY_FLOAT64:
        _apply_cals_typed<double>(array, cals);
        break;
    default:
        throw g3supertimestream_exception("Unexpected dtype!");
    }
}

void G3SuperTimestream::Calibrate(std::vector<double> rescale)
{
    if (rescale.size() != names.size())
        throw g3supertimestream_exception(
            "Rescale vector has unexpected length.");
    if (float_mode) {
        // Modification to the calibration.
        if (array)
            _apply_cals(array, rescale);
        for (size_t i=0; i<quanta.size(); i++)
            quanta[i] *= rescale[i];
    } else {
        // Transition to float_mode.  If holding integer
        // array, convert it.
        if (dataful) {
            switch(desc.type_num) {
            case NPY_INT32:
                desc.type_num = NPY_FLOAT32;
                break;
            case NPY_INT64:
                desc.type_num = NPY_FLOAT64;
                break;
            default:
                throw g3supertimestream_exception("Unexpected dtype!");
            }
            if (array) {
                auto new_array = (PyArrayObject*)PyArray_FromAny(
                    (PyObject*)array, PyArray_DescrFromType(desc.type_num),
                    0, 0, NPY_ARRAY_FORCECAST, NULL);
                if (new_array == nullptr)
                    throw g3supertimestream_exception(
                        "Failed to allocate float array.");
                _apply_cals(new_array, rescale);
                Py_DECREF(array);
                array = new_array;
            }
        }
        float_mode = true;
        quanta = rescale;
    }
}

G3SuperTimestream::G3SuperTimestream()
{
    array = nullptr;
    ablob = nullptr;
    float_mode = false;
    dataful = false;
    if (!PyArray_EquivByteorders(NPY_NATIVE, NPY_LITTLE))
        throw g3supertimestream_exception(
            "This class hasn't been trained on BIG-endian machines.!");

    // Set defaults.
    Options(1, 1, 30, -1, -1);
}

G3SuperTimestream::~G3SuperTimestream()
{
    if (array != nullptr)
        Py_XDECREF(array);
    if (ablob != nullptr) {
        delete ablob->buf;
        delete ablob;
    }
}

G3SuperTimestream::G3SuperTimestream(
    const G3VectorString &names_, const G3VectorTime &times_)
    : G3SuperTimestream()
{
    names = G3VectorString(names_);
    times = G3VectorTime(times_);
}

static
void safe_set_data(G3SuperTimestream &self, const bp::object object_in);
static
void safe_set_quanta(G3SuperTimestream &self, std::vector<double> quanta);

G3SuperTimestream::G3SuperTimestream(
    const G3VectorString &names_, const G3VectorTime &times_,
    const bp::object &data_)
    : G3SuperTimestream(names_, times_)
{
    safe_set_data(*this, data_);
}

G3SuperTimestream::G3SuperTimestream(
    const G3VectorString &names_, const G3VectorTime &times_,
    const bp::object &data_, const std::vector<double> &quanta_)
    : G3SuperTimestream(names_, times_)
{
    safe_set_quanta(*this, quanta_);
    safe_set_data(*this, data_);
}

static
void safe_set_times(G3SuperTimestream &self, G3VectorTime _times)
{
    // Only allow this if it doesn't upset consistency.  We will
    // assume that, coming in, we're internally consistent.
    if (_times.size() != self.times.size() && self.times.size() != 0) {
        std::ostringstream s;
        s << "Cannot set .times because it conflicts with "
          << "the established number of samples (" << self.times.size()
          << ").";
        throw g3supertimestream_exception(s.str());
    }
    self.times = _times;
}

static
void safe_set_names(G3SuperTimestream &self, G3VectorString _names)
{
    // Only allow this if it doesn't upset consistency.  We will
    // assume that, coming in, we're internally consistent.
    if (_names.size() != self.names.size() && self.names.size() != 0) {
        std::ostringstream s;
        s << "Cannot set .names because it conflicts with "
          << "the established number of channels (" << self.names.size()
          << ").";
        throw g3supertimestream_exception(s.str());
    }
    self.names = _names;
}

static
void safe_set_data(G3SuperTimestream &self, const bp::object object_in)
{
    // Note this function, as invoked here, might return a
    // reference or create a new array.
    PyObject *ob = PyArray_FromAny(object_in.ptr(), NULL, 0, 0, 0, NULL);
    if (ob == NULL)
        throw g3supertimestream_exception("Could not decode array.");

    PyArrayObject *_array = reinterpret_cast<PyArrayObject*>(ob);

    if (PyArray_NDIM(_array) != 2) {
        Py_XDECREF(ob);
        throw g3supertimestream_exception("Bad ndim.");
    }
    if ((size_t)PyArray_DIMS(_array)[0] != self.names.size()) {
        Py_XDECREF(ob);
        throw g3supertimestream_exception("Bad shape[0].");
    }
    if ((size_t)PyArray_DIMS(_array)[1] != self.times.size()) {
        Py_XDECREF(ob);
        throw g3supertimestream_exception("Bad shape[1].");
    }
    if (!PyArray_EquivByteorders(PyArray_DESCR(_array)->byteorder, NPY_LITTLE)) {
        //There are other ways to deal with endianness.
        Py_XDECREF(ob);
        throw g3supertimestream_exception("Bad endianness.");
    }

    bool is_floaty = false;
    switch(PyArray_TYPE(_array)) {
    case NPY_FLOAT32:
    case NPY_FLOAT64:
        is_floaty = true;
        break;
    case NPY_INT32:
    case NPY_INT64:
        break;
    default:
        Py_XDECREF(ob);
        throw g3supertimestream_exception("Forbidden dtype.");
    }

    if (is_floaty) {
        // quanta has to be set already.
        if (self.quanta.size() != (size_t)PyArray_DIMS(_array)[0])
            throw g3supertimestream_exception(
                "User must set .quanta before loading float array.");
    } else {
        if (self.quanta.size() != 0)
            throw g3supertimestream_exception(
                "The .quanta must be empty when loading integer array.");
    }

    // Clear cached array or compressed data.
    if (self.array) {
        Py_XDECREF((PyObject*)self.array);
        self.array = nullptr;
    }
    if (self.ablob) {
        delete self.ablob->buf;
        delete self.ablob;
        self.ablob = nullptr;
    }

    self.dataful = true;
    self.float_mode = is_floaty;

    self.desc.ndim = PyArray_NDIM(_array);
    self.desc.type_num = PyArray_TYPE(_array);

    self.desc.nbytes = PyArray_NBYTES(_array);
    for (int i=0; i<self.desc.ndim; i++)
        self.desc.shape[i] = PyArray_DIMS(_array)[i];

    self.array = _array;
}

static
bp::object safe_get_data(G3SuperTimestream &self)
{
    if (self.array == nullptr)
        self.Decode();
    return bp::object(bp::handle<>(bp::borrowed(reinterpret_cast<PyObject*>(self.array))));
}

static
bp::object safe_get_dtype(G3SuperTimestream &self)
{
    PyObject *d;
    if (self.array == nullptr) {
        d = reinterpret_cast<PyObject*>(
            PyArray_DescrFromType(self.desc.type_num));
    } else {
        d = reinterpret_cast<PyObject*>(PyArray_DESCR(self.array));
        Py_XINCREF(d);
    }
    return bp::object(bp::handle<>(d));
}

static
bp::object safe_get_quanta(G3SuperTimestream &self)
{
    if (!self.float_mode)
        return bp::object();

    npy_intp shape[1] = {(npy_intp)self.quanta.size()};
    auto output = (PyArrayObject *)PyArray_SimpleNew(1, shape, NPY_FLOAT64);
    memcpy(PyArray_DATA(output), &self.quanta[0], shape[0] * sizeof(self.quanta[0]));
    return bp::object(bp::handle<>(reinterpret_cast<PyObject*>(output)));
}

static
void safe_set_quanta(G3SuperTimestream &self, std::vector<double> quanta)
{
    // Only allowed to set quanta directly if data is not present.
    if (!self.dataful)
        self.Calibrate(quanta);
    else
        throw g3supertimestream_exception(
            "The .quanta cannot be set directly once .data is set.  Use .calibrate().");
}

//Copies data out of a flat buffer, creating the numpy array along the way.
bool G3SuperTimestream::SetDataFromBuffer(void* buf, int ndim, int shape[], int typenum,
                     std::pair<int,int> sample_range)
{
    if (ndim != 2)
        throw g3supertimestream_exception(
            "2d arrays only please");

    // Create a new numpy array for this, allowing for slice in
    // second dimension..
    int n_samps = sample_range.second - sample_range.first;
    npy_intp shape_[2] = {shape[0], n_samps};
    auto array_ = (PyArrayObject*)PyArray_EMPTY(ndim, shape_, typenum, 0);
    bp::object array_ob =
        bp::object(bp::handle<>((reinterpret_cast<PyObject*>(array_))));

    for (int i=0; i<shape[0]; i++) {
        memcpy((char*)PyArray_DATA(array_) + PyArray_STRIDES(array_)[0] * i,
               (char*)buf + (sample_range.first + shape[1] * i) * PyArray_ITEMSIZE(array_),
               PyArray_STRIDES(array_)[1] * n_samps);
    }

    safe_set_data(*this, array_ob);
    return true;
}


// Assist with testing the pure C++ interface
static
G3SuperTimestreamPtr test_cxx_interface(int nsamps, int first, int second)
{
    int shape[2] = {3, nsamps};
    int typenum = NPY_INT32;
    std::vector<int32_t> buf(shape[0] * shape[1]);
    std::fill(buf.begin(), buf.end(), 0);

    auto times = G3VectorTime();
    for (int i=first; i<second; i++)
        times.push_back(G3Time::Now());
    const char *chans_[] = {"a", "b", "c"};
    auto chans = G3VectorString(chans_, std::end(chans_));

    auto ts = G3SuperTimestreamPtr(
        new G3SuperTimestream(chans, times));

    for (int i=first; i<second; i++)
        buf[i] = 77;

    ts->SetDataFromBuffer((void*)buf.data(), 2, shape, typenum,
                  std::pair<int,int>(first, second));

    return ts;
}


G3_SPLIT_SERIALIZABLE_CODE(G3SuperTimestream);

static void translate_ValueError(g3supertimestream_exception const& e)
{
    PyErr_SetString(PyExc_ValueError, e.msg_for_python().c_str());
}

static void* _so_import_array() {
    import_array();
    return NULL;
}

PYBINDINGS("so")
{
    using namespace boost::python;

    bp::docstring_options local_docstring_options(true, false);

    // ndarray import
    _so_import_array();

    EXPORT_FRAMEOBJECT(G3SuperTimestream, init<>(), "G3SuperTimestream()\n\n"
        "  Construct with names and times uninitialized.")
        .def(bp::init<const G3VectorString&, const G3VectorTime&>(
                 "G3SuperTimestream(names, times)\n"
                 "  Construct with names and times initialized."))
        .def(bp::init<const G3VectorString&, const G3VectorTime&,
             bp::object&>(
                 "G3SuperTimestream(names, times, data)\n\n"
                 "  Construct with integer data."))
        .def(bp::init<const G3VectorString&, const G3VectorTime&,
             bp::object&, const std::vector<double>&>(
                 "G3SuperTimestream(names, times, data, quanta)\n\n"
                 "  Construct with float data."))
        .add_property("times", &G3SuperTimestream::times, &safe_set_times,
                  "Vector of timestamps (G3VectorTime).")
        .add_property("names", &G3SuperTimestream::names, &safe_set_names,
                  "Vector of channel names.")
        .add_property("data", &safe_get_data, &safe_set_data,
                  "Data array.")
        .add_property("quanta", &safe_get_quanta, &safe_set_quanta,
                  "Quanta (if float mode).")
        .add_property("dtype", &safe_get_dtype, "Numpy dtype of enclosed array.")
        .def("encode", &G3SuperTimestream::Encode, "Compress.")
        .def("decode", &G3SuperTimestream::Decode, "Decompress.")
        .def("extract", &G3SuperTimestream::Extract,
             (bp::arg("dest"), bp::arg("dest_indices")=bp::object(),
              bp::arg("src_indices")=bp::object(),
              bp::arg("start")=0, bp::arg("stop")=-1),
             "Decompress data subset into an array.")
        .def("calibrate", &G3SuperTimestream::Calibrate,
             "calibrate(cal_factors)\n\n"
                     "Apply per-channel scale factors.  Note this puts you in float mode\n"
             "(if you're not already) and modifies both .quanta and .data.")
        .def("options", &G3SuperTimestream::Options,
             (bp::arg("enable")=-1,
              bp::arg("flac_level")=-1, bp::arg("bz2_workFactor")=-1,
              bp::arg("data_algo")=-1, bp::arg("times_algo")=-1),
             "Set compression options.")
        ;
    register_pointer_conversions<G3SuperTimestream>();

    bp::register_exception_translator<g3supertimestream_exception>(&translate_ValueError);
    bp::def("test_g3super", test_cxx_interface);
}
