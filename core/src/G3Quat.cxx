#include <pybindings.h>
#include <container_pybindings.h>
#include <G3Quat.h>
#include <G3Map.h>

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


G3_SERIALIZABLE_CODE(G3VectorQuat);

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

template <>
G3VectorQuatPtr
container_from_object(boost::python::object v)
{
        G3VectorQuatPtr x(new (G3VectorQuat));
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
	     .def(self /= self)
	     .def(self /= double())
	     .def("__str__", quat_str)
	     .def("__repr__", quat_repr)
	     .def("dot3", dot3, "Dot product of last three entries")
	     .def("cross3", cross3, "Cross product of last three entries")
	;
	register_vector_of<quat>("QuatVector");
	object vq =
	    register_g3vector<quat>("G3VectorQuat", "List of quaternions")
	     .def(self * double())
	     .def(double() * self)
	     .def(self * self)
	     .def(self * quat())
	     .def(quat() * self)
	     .def(self *= double())
	     .def(self *= quat())
	     .def(self *= self)
	     .def(self / double())
	     .def(self /= double())
	     .def(self / self)
	     .def(self /= self)
	     .def(pow(self, double()))
	     .def(pow(self, int()));
	PyTypeObject *vqclass = (PyTypeObject *)vq.ptr();
	vectorquat_bufferprocs.bf_getbuffer = G3VectorQuat_getbuffer;
	vqclass->tp_as_buffer = &vectorquat_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vqclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

}
