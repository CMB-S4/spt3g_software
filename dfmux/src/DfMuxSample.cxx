#include <pybindings.h>
#include <container_pybindings.h>

#include <dfmux/DfMuxSample.h>

template <class A> void DfMuxSample::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("samples", base_class<std::vector<int32_t> >(this));
	ar & make_nvp("time", Timestamp);
}

PYBINDINGS("dfmux", scope)
{
	register_g3vector<DfMuxSample>(scope, "DfMuxSample",
	  "Samples from all channels on one readout module, stored with I and "
	  "Q interleaved, such that the first element is channel 1 I, followed "
	  "by channel 1 Q, followed by channel 2 I, etc.")
	    .def(py::init<G3TimeStamp, int>(), py::arg("time"), py::arg("nsamples"))
	    .def_readwrite("Timestamp", &DfMuxSample::Timestamp)
	;
}

G3_SERIALIZABLE_CODE(DfMuxSample);

