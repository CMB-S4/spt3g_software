#include <pybindings.h>
#include <container_pybindings.h>
#include <G3Quat.h>
#include <G3Map.h>
#include <G3Units.h>

// Quaternion utilities

quat
cross3(quat u, quat v)
{
	// Computes Euclidean cross product from the last three entries in the
	// quaternion
	return quat( 
	    0, 
	    u.R_component_3()*v.R_component_4() - (u.R_component_4()*v.R_component_3()),
	    u.R_component_4()*v.R_component_2() - (u.R_component_2()*v.R_component_4()),
	    u.R_component_2()*v.R_component_3() - (u.R_component_3()*v.R_component_2()));
}

double
dot3(quat a, quat b)
{
	// Computes Euclidean dot product from the last three entries in the
	// quaternion
	return (a.R_component_2()*b.R_component_2() +
		a.R_component_3()*b.R_component_3() +
		a.R_component_4()*b.R_component_4());
}

static double
_abs(const quat &a)
{
	return sqrt(norm(a));
}

static G3VectorDouble
_vabs(const G3VectorQuat &a)
{
	G3VectorDouble out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
                out[i] = _abs(a[i]);
	return out;
}

namespace boost {
namespace math {
quat
operator ~(quat a)
{
	return conj(a);
}
};
};

G3VectorQuat
operator ~(const G3VectorQuat &a)
{
	G3VectorQuat out(a.size());
        for (unsigned i = 0; i < a.size(); i++)
		out[i] = conj(a[i]);
	return out;
}

G3VectorQuat
operator *(const G3VectorQuat &a, double b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3VectorQuat
operator *(double b, const G3VectorQuat &a)
{
	return a*b;
}

G3VectorQuat &
operator *=(G3VectorQuat &a, double b)
{
	for (quat &i: a)
		i *= b;
	return a;
}

G3VectorQuat
operator /(const G3VectorQuat &a, double b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3VectorQuat
operator /(double a, const G3VectorQuat &b)
{
	G3VectorQuat out(b.size());
	for (unsigned i = 0; i < b.size(); i++)
		out[i] = a/b[i];
	return out;
}

G3VectorQuat
operator /(const G3VectorQuat &a, const quat &b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3VectorQuat
operator /(const quat &a, const G3VectorQuat &b)
{
	G3VectorQuat out(b.size());
	for (unsigned i = 0; i < b.size(); i++)
		out[i] = a/b[i];
	return out;
}

G3VectorQuat
operator /(const G3VectorQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b[i];
	return out;
}

G3VectorQuat &
operator /=(G3VectorQuat &a, double b)
{
	for (quat &i: a)
		i /= b;
	return a;
}

G3VectorQuat &
operator /=(G3VectorQuat &a, const quat &b)
{
	for (unsigned i = 0; i < a.size(); i++)
		a[i] /= b;
	return a;
}

G3VectorQuat &
operator /=(G3VectorQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	for (unsigned i = 0; i < a.size(); i++)
		a[i] /= b[i];
	return a;
}

G3VectorQuat
operator *(const G3VectorQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b[i];
	return out;
}

G3VectorQuat &
operator *=(G3VectorQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	for (unsigned i = 0; i < a.size(); i++)
		a[i] *= b[i];
	return a;
}

G3VectorQuat
operator *(const G3VectorQuat &a, quat b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3VectorQuat
operator *(quat b, const G3VectorQuat &a)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = b*a[i];
	return out;
}

G3VectorQuat &
operator *=(G3VectorQuat &a, quat b)
{
	for (quat &i: a)
		i *= b;
	return a;
}

G3VectorQuat
pow(const G3VectorQuat &a, double b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = pow(a[i], b);
	return out;
}

G3VectorQuat
pow(const G3VectorQuat &a, int b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = pow(a[i], b);
	return out;
}

double G3TimestreamQuat::GetSampleRate() const
{
	return double(size() - 1)/double(stop.time - start.time);
}

std::string
G3TimestreamQuat::Description() const
{
	std::ostringstream desc;
	desc.precision(1);
	desc << std::fixed;
	desc << size() << " quaternions at " << GetSampleRate()/G3Units::Hz <<
	    " Hz";

	return desc.str();
}

template <class A>
void G3TimestreamQuat::serialize(A &ar, const unsigned v)
{
        G3_CHECK_VERSION(v);

        ar & cereal::make_nvp("G3VectorQuat",
            cereal::base_class<G3VectorQuat>(this));
        ar & cereal::make_nvp("start", start);
        ar & cereal::make_nvp("stop", stop);
}

// Following duplicates above precisely because of return types. This should
// probably use templates or pre-processor macros but that seemed to
// take almost as much work as copying, especially given the need to copy
// start/stop times in the non-in-place operators.
G3TimestreamQuat
operator ~(const G3TimestreamQuat &a)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
        for (unsigned i = 0; i < a.size(); i++)
		out[i] = conj(a[i]);
	return out;
}

G3TimestreamQuat
operator *(const G3TimestreamQuat &a, double b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3TimestreamQuat
operator *(double b, const G3TimestreamQuat &a)
{
	return a*b;
}

G3TimestreamQuat &
operator *=(G3TimestreamQuat &a, double b)
{
	for (quat &i: a)
		i *= b;
	return a;
}

G3TimestreamQuat
operator /(const G3TimestreamQuat &a, double b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3TimestreamQuat
operator /(double a, const G3TimestreamQuat &b)
{
	G3TimestreamQuat out(b.size());
	out.start = b.start; out.stop = b.stop;
	for (unsigned i = 0; i < b.size(); i++)
		out[i] = a/b[i];
	return out;
}

G3TimestreamQuat
operator /(const G3TimestreamQuat &a, const quat &b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3TimestreamQuat
operator /(const quat &a, const G3TimestreamQuat &b)
{
	G3TimestreamQuat out(b.size());
	out.start = b.start; out.stop = b.stop;
	for (unsigned i = 0; i < b.size(); i++)
		out[i] = a/b[i];
	return out;
}

G3TimestreamQuat
operator /(const G3TimestreamQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b[i];
	return out;
}

G3TimestreamQuat &
operator /=(G3TimestreamQuat &a, double b)
{
	for (quat &i: a)
		i /= b;
	return a;
}

G3TimestreamQuat &
operator /=(G3TimestreamQuat &a, const quat &b)
{
	for (unsigned i = 0; i < a.size(); i++)
		a[i] /= b;
	return a;
}

G3TimestreamQuat &
operator /=(G3TimestreamQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	for (unsigned i = 0; i < a.size(); i++)
		a[i] /= b[i];
	return a;
}

G3TimestreamQuat
operator *(const G3TimestreamQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b[i];
	return out;
}

G3TimestreamQuat &
operator *=(G3TimestreamQuat &a, const G3VectorQuat &b)
{
	g3_assert(a.size() == b.size());
	for (unsigned i = 0; i < a.size(); i++)
		a[i] *= b[i];
	return a;
}

G3TimestreamQuat
operator *(const G3TimestreamQuat &a, quat b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3TimestreamQuat
operator *(quat b, const G3TimestreamQuat &a)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = b*a[i];
	return out;
}

G3TimestreamQuat &
operator *=(G3TimestreamQuat &a, quat b)
{
	for (quat &i: a)
		i *= b;
	return a;
}

G3TimestreamQuat
pow(const G3TimestreamQuat &a, double b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = pow(a[i], b);
	return out;
}

G3TimestreamQuat
pow(const G3TimestreamQuat &a, int b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = pow(a[i], b);
	return out;
}


G3_SERIALIZABLE_CODE(G3VectorQuat);
G3_SERIALIZABLE_CODE(G3TimestreamQuat);

namespace {
static int
G3VectorQuat_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	bp::handle<> self(bp::borrowed(obj));
	bp::object selfobj(self);
	G3VectorQuatPtr q = bp::extract<G3VectorQuatPtr>(selfobj)();

	view->obj = obj;
	view->buf = (void*)&(*q)[0];
	view->len = q->size() * sizeof(double) * 4;
	view->readonly = 0;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;

	// XXX: following leaks small amounts of memory!
	view->shape = new Py_ssize_t[2];
	view->strides = new Py_ssize_t[2];

	view->ndim = 2;
	view->shape[0] = q->size();
	view->shape[1] = 4;
	view->strides[0] = view->shape[1]*view->itemsize;
	view->strides[1] = view->itemsize;

	view->suboffsets = NULL;

	Py_INCREF(obj);

	return 0;
}

static PyBufferProcs vectorquat_bufferprocs;
static PyBufferProcs timestreamquat_bufferprocs;

static std::string
quat_str(const quat &q)
{
	std::ostringstream oss;
	oss << q;
	return oss.str();
}

static std::string
quat_repr(const quat &q)
{
	std::ostringstream oss;
	oss << "spt3g.core.quat" << q;
	return oss.str();
}
}

template <typename T>
boost::shared_ptr<T>
quat_vec_container_from_object(boost::python::object v)
{
        boost::shared_ptr<T> x(new (T));
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (view.ndim != 2 || view.shape[1] != 4) {
			boost::python::container_utils::extend_container(*x, v);
		} else if (strcmp(view.format, "d") == 0) {
			x->resize(view.shape[0]);
			memcpy((void *)&(*x)[0], view.buf, view.len);
		} else if (strcmp(view.format, "f") == 0) {
			x->resize(view.shape[0]);
			for (size_t i = 0; i < view.shape[0]; i++)
				(*x)[i] = quat(
				    ((float *)view.buf)[4*i + 0],
				    ((float *)view.buf)[4*i + 1],
				    ((float *)view.buf)[4*i + 2],
				    ((float *)view.buf)[4*i + 3]);
		} else if (strcmp(view.format, "i") == 0) {
			x->resize(view.shape[0]);
			for (size_t i = 0; i < view.shape[0]; i++)
				(*x)[i] = quat(
				    ((int *)view.buf)[4*i + 0],
				    ((int *)view.buf)[4*i + 1],
				    ((int *)view.buf)[4*i + 2],
				    ((int *)view.buf)[4*i + 3]);
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

	return x;
}

template <>
G3VectorQuatPtr container_from_object(boost::python::object v)
{
	return quat_vec_container_from_object<G3VectorQuat>(v);
}

template <>
G3TimestreamQuatPtr container_from_object(boost::python::object v)
{
	return quat_vec_container_from_object<G3TimestreamQuat>(v);
}

static int
G3TimestreamQuat_nsamples(const G3TimestreamQuat &r)
{
        return r.size();
}

PYBINDINGS("core")
{
	using namespace boost::python;

	class_<quat>("quat",
	    "Representation of a quaternion. Data in a,b,c,d.",
	    init<double, double, double, double>())
	     .add_property("a", &quat::R_component_1)
	     .add_property("b", &quat::R_component_2)
	     .add_property("c", &quat::R_component_3)
	     .add_property("d", &quat::R_component_4)
	     .def(~self)
	     .def(self == self)
	     .def(self != self)
	     .def(self + self)
	     .def(self += self)
	     .def(self - self)
	     .def(self -= self)
	     .def(self * self)
	     .def(self * double())
	     .def(double() * self)
	     .def(self *= self)
	     .def(self *= double())
	     .def(pow(self, double()))
	     .def(pow(self, long()))
	     .def(self / self)
	     .def(self / double())
	     .def(double() / self)
	     .def(self /= self)
	     .def(self /= double())
	     .def("__abs__", _abs)
	     .def("__str__", quat_str)
	     .def("__repr__", quat_repr)
	     .def("dot3", dot3, "Dot product of last three entries")
	     .def("cross3", cross3, "Cross product of last three entries")
	;
	register_vector_of<quat>("QuatVector");
	object vq =
	    register_g3vector<quat>("G3VectorQuat",
	     "List of quaternions. Convertible to a 4xN numpy array. "
	     "Arithmetic operations on this object are fast and provide "
	     "results given proper quaternion math rather than "
	     "element-by-element numpy-ish results.")
	     .def(~self)
	     .def(self * double())
	     .def(double() * self)
	     .def(self * self)
	     .def(self * quat())
	     .def(quat() * self)
	     .def(self *= double())
	     .def(self *= quat())
	     .def(self *= self)
	     .def(self / double())
	     .def(double() / self)
	     .def(self /= double())
	     .def(self / self)
	     .def(self /= self)
	     .def(self / quat())
	     .def(self /= quat())
	     .def(quat() / self)
	     .def(pow(self, double()))
	     .def(pow(self, int()))
	     .def("__abs__", _vabs);
	PyTypeObject *vqclass = (PyTypeObject *)vq.ptr();
	vectorquat_bufferprocs.bf_getbuffer = G3VectorQuat_getbuffer;
	vqclass->tp_as_buffer = &vectorquat_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vqclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	object tq =
	    class_<G3TimestreamQuat, bases<G3VectorQuat>, G3TimestreamQuatPtr>(
	      "G3TimestreamQuat",
	      "Timestream of quaternions. Identical to a G3VectorQuat except "
	      "for the addition of start and stop times.",
	      init<>()
             )
	     .def("__init__", make_constructor(container_from_object<G3TimestreamQuat>))
	     .def(boost::python::init<const G3TimestreamQuat &>())
	     .def_pickle(g3frameobject_picklesuite<G3TimestreamQuat>())
	     .def(~self)
	     .def(self * double())
	     .def(double() * self)
	     .def(self * G3VectorQuat())
	     .def(self * quat())
	     .def(quat() * self)
	     .def(self *= double())
	     .def(self *= quat())
	     .def(self *= G3VectorQuat())
	     .def(self / double())
	     .def(double() / self)
	     .def(self /= double())
	     .def(self / G3VectorQuat())
	     .def(self /= G3VectorQuat())
	     .def(self / quat())
	     .def(self /= quat())
	     .def(quat() / self)
	     .def(pow(self, double()))
	     .def(pow(self, int()))
	     .def("__abs__", _vabs)
	    .def_readwrite("start", &G3TimestreamQuat::start,
	      "Time of the first sample in the time stream")
	    .def_readwrite("stop", &G3TimestreamQuat::stop,
	      "Time of the final sample in the timestream")
	    .add_property("sample_rate", &G3TimestreamQuat::GetSampleRate,
	       "Computed sample rate of the timestream.")
	    .add_property("n_samples", &G3TimestreamQuat_nsamples,
	      "Number of samples in the timestream. Equivalent to len(ts)")
	;
	PyTypeObject *tqclass = (PyTypeObject *)tq.ptr();
	timestreamquat_bufferprocs.bf_getbuffer = G3VectorQuat_getbuffer;
	tqclass->tp_as_buffer = &timestreamquat_bufferprocs;
#if PY_MAJOR_VERSION < 3
	tqclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	scitbx::boost_python::container_conversions::from_python_sequence<G3TimestreamQuat, scitbx::boost_python::container_conversions::variable_capacity_policy>();
	register_pointer_conversions<G3TimestreamQuat>();
	implicitly_convertible<G3TimestreamQuatPtr, G3VectorQuatPtr>();
	implicitly_convertible<G3TimestreamQuatPtr, G3VectorQuatConstPtr>();
}
