#include <pybindings.h>
#include <G3Module.h>


class G3InfiniteSource : public G3Module {
public:
	G3InfiniteSource(G3Frame::FrameType type = G3Frame::None, int n = -1) : type_(type), n_(n), sofar_(0) {}
	
	void Process(G3FramePtr frame, std::deque<G3FramePtr> &out) {
		sofar_++;
		if (n_ < 0 || sofar_ <= n_)
			out.push_back(G3FramePtr(new G3Frame(type_)));
	}
private:
	G3Frame::FrameType type_;
	int n_, sofar_;
};

PYBINDINGS("core", scope) {
	register_g3module<G3InfiniteSource>(scope, "G3InfiniteSource",
	    "Emits infinite frames, up to an optional maximum number n")
	    .def(py::init<>())
	    .def(py::init<G3Frame::FrameType, int>(),
	        py::arg("type")=G3Frame::None, py::arg("n")=-1);
};
