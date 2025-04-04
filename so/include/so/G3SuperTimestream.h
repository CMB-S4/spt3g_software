#pragma once

#include <so/so3g_numpy.h>
#include <so/numpy_assist.h>
#include <G3Frame.h>
#include <G3Map.h>

#include <exception>
#include <stdint.h>

using namespace std;


class G3SuperTimestream : public G3FrameObject {
    // Storage for a 2d array with shape (n_dets, n_samps) along
    // with the vector of associated timestamps (n_samps), and the
    // vector of channel names (n_dets).  Serializes with lossless
    // compression of int32 and int64, and controllable
    // quantization and compression of float32 and lfoat64.
public:
    G3SuperTimestream();
    ~G3SuperTimestream();

    G3SuperTimestream(const G3VectorString &names_, const G3VectorTime &times_);
    G3SuperTimestream(const G3VectorString &names_, const G3VectorTime &times_,
              const bp::object &data);
    G3SuperTimestream(const G3VectorString &names_, const G3VectorTime &times_,
              const bp::object &data_, const std::vector<double> &quanta_);

    // This object contains pointers to memory that it
    // allocated... and they're freed on destruction.  The
    // responsible thing to do in such circumstances is to delete
    // or modify the move and copy constructors.  But if we do
    // that, boost python complains about something.  In G3 you
    // avoid move/copy by wrapping all instances in a G3xxxPtr,
    // and essentially pass the object around by reference, so
    // that's what we'll do.  But beware, there's nothing
    // preventing you from returning G3SuperTimestream directly
    // from a function (except the inevitable segfault).

    string Description() const;
    string Summary() const;

    bool Extract(bp::object dest, bp::object dest_indices, bp::object src_indices,
             int start, int stop);
    bool Encode();
    bool Decode();
    void Calibrate(vector<double> rescale);
    int Options(int enable=-1,
            int flac_level=-1, int bz2_workFactor=-1,
                    int data_algo=-1, int times_algo=-1);

    // Interface for C++...
    bool SetDataFromBuffer(void* buf, int ndim, int shape[], int typenum,
                   std::pair<int,int> sample_range);


    template <class A> void load(A &ar, unsigned v);
    template <class A> void save(A &ar, unsigned v) const;

    struct array_desc {
        npy_intp type_num;
        npy_intp ndim;
        npy_intp shape[32];
        npy_intp nbytes;
    };

    // Container for the compressed data.
    struct array_blob {
        int size;
        char *buf;
        int count;
        vector<int> offsets;
    };

    enum algos {
        ALGO_NONE = 0,
        ALGO_DO_FLAC = (1 << 0),
        ALGO_DO_BZ = (1 << 1),
        ALGO_DO_CONST = (1 << 2)
    };

    struct options_type {
        int8_t times_algo;
        int8_t data_algo;
        int8_t flac_level;
        int8_t bz2_workFactor;
    } options;

    G3VectorTime times;
    G3VectorString names;

    bool float_mode;
    bool dataful;
    vector<double> quanta;
    struct array_desc desc;

    PyArrayObject *array;
    struct array_blob *ablob;
};

// This specialization tells cereal to use G3SuperTimestream::load/save
// and not the base class' load/save.
namespace cereal {
    template <class A> struct specialize<
        A, G3SuperTimestream, cereal::specialization::member_load_save> {};
}

G3_POINTERS(G3SuperTimestream);
G3_SERIALIZABLE(G3SuperTimestream, 0);

class g3supertimestream_exception : std::exception
{
    // Exception raised when internal validity checks fail.  This will
    // also be mapped to some particular Python exception type.
public:
    std::string text;
    g3supertimestream_exception(std::string text) :	text{text} {};

    std::string msg_for_python() const throw() {
        return text;
    }
};
