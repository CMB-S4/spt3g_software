#include <pybindings.h>
#include <serialization.h>

#include <dfmux/DfMuxSample.h>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>

template <class A> void DfMuxSample::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("samples", base_class<std::vector<int32_t> >(this));
	ar & make_nvp("time", Timestamp);
}

PYBINDINGS("dfmux")
{
	namespace bp = boost::python;

	bp::class_<DfMuxSample, bp::bases<G3FrameObject, std::vector<int32_t> >,
	  DfMuxSamplePtr, boost::noncopyable>("DfMuxSample",
	  "Samples from all channels on one readout module, stored with I and "
	  "Q interleaved, such that the first element is channel 1 I, followed "
	  "by channel 1 Q, followed by channel 2 I, etc.",
	  bp::init<G3TimeStamp, int>(bp::args("time", "nsamples")))
	    .def_readwrite("Timestamp", &DfMuxSample::Timestamp)
	    .def_pickle(g3frameobject_picklesuite<DfMuxSample>())
	;
	register_pointer_conversions<DfMuxSample>();
}

G3_SERIALIZABLE_CODE(DfMuxSample);

