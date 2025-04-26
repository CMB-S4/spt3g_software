#include <pybindings.h>
#include <container_pybindings.h>
#include <G3Quat.h>
#include <G3Map.h>
#include <G3Units.h>
#include <G3Timestream.h>

// Quaternion utilities

std::ostream&
operator<<(std::ostream& os, const Quat &q)
{
	os << "(" << q.a() << ", " << q.b() << ", " << q.c() << ", " << q.d() << ")";
	return os;
}

template <class A>
void
Quat::serialize(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("a", a_);
	ar & cereal::make_nvp("b", b_);
	ar & cereal::make_nvp("c", c_);
	ar & cereal::make_nvp("d", d_);
}

std::string
G3Quat::Description() const
{
	std::ostringstream oss;
	oss << value;
	return oss.str();
}

template <class A>
void
G3Quat::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("value", value);
}

double
Quat::real() const
{
	return a_;
}

Quat
Quat::unreal() const
{
	if (!a_)
		return *this;
	return Quat(0, b_, c_, d_);
}

Quat
Quat::conj() const
{
	return Quat(a_, -b_, -c_, -d_);
}

double
Quat::norm() const
{
	return a_ * a_ + b_ * b_ + c_ * c_ + d_ * d_;
}

double
Quat::vnorm() const
{
	return b_ * b_ + c_ * c_ + d_ * d_;
}

double
Quat::abs() const
{
	return sqrt(norm());
}

Quat
Quat::operator -() const
{
	return Quat(-a_, -b_, -c_, -d_);
}

Quat
Quat::operator ~() const
{
	return conj();
}

Quat &
Quat::operator +=(const Quat &rhs)
{
	a_ += rhs.a_;
	b_ += rhs.b_;
	c_ += rhs.c_;
	d_ += rhs.d_;
	return *this;
}

Quat &
Quat::operator -=(const Quat &rhs)
{
	a_ -= rhs.a_;
	b_ -= rhs.b_;
	c_ -= rhs.c_;
	d_ -= rhs.d_;
	return *this;
}

Quat &
Quat::operator *=(double rhs)
{
	a_ *= rhs;
	b_ *= rhs;
	c_ *= rhs;
	d_ *= rhs;
	return *this;
}

Quat &
Quat::operator *=(const Quat &rhs)
{
	double a = a_ * rhs.a_ - b_ * rhs.b_ - c_ * rhs.c_ - d_ * rhs.d_;
	double b = a_ * rhs.b_ + b_ * rhs.a_ + c_ * rhs.d_ - d_ * rhs.c_;
	double c = a_ * rhs.c_ - b_ * rhs.d_ + c_ * rhs.a_ + d_ * rhs.b_;
	double d = a_ * rhs.d_ + b_ * rhs.c_ - c_ * rhs.b_ + d_ * rhs.a_;
	a_ = a;
	b_ = b;
	c_ = c;
	d_ = d;
	return *this;
}

Quat &
Quat::operator /=(double rhs)
{
	a_ /= rhs;
	b_ /= rhs;
	c_ /= rhs;
	d_ /= rhs;
	return *this;
}

Quat &
Quat::operator /=(const Quat &rhs)
{
	double n = rhs.norm();
	double a =  a_ * rhs.a_ + b_ * rhs.b_ + c_ * rhs.c_ + d_ * rhs.d_;
	double b = -a_ * rhs.b_ + b_ * rhs.a_ - c_ * rhs.d_ + d_ * rhs.c_;
	double c = -a_ * rhs.c_ + b_ * rhs.d_ + c_ * rhs.a_ - d_ * rhs.b_;
	double d = -a_ * rhs.d_ - b_ * rhs.c_ + c_ * rhs.b_ + d_ * rhs.a_;
	a_ = a / n;
	b_ = b / n;
	c_ = c / n;
	d_ = d / n;
	return *this;
}

Quat
Quat::operator +(const Quat &rhs) const
{
	return Quat(a_ + rhs.a_, b_ + rhs.b_, c_ + rhs.c_, d_ + rhs.d_);
}

Quat
Quat::operator -(const Quat &rhs) const
{
	return Quat(a_ - rhs.a_, b_ - rhs.b_, c_ - rhs.c_, d_ - rhs.d_);
}

Quat
Quat::operator *(double rhs) const
{
	return Quat(a_ * rhs, b_ * rhs, c_ * rhs, d_ * rhs);
}

Quat
Quat::operator *(const Quat &rhs) const
{
	Quat out(*this);
	out *= rhs;
	return out;
}

Quat
operator *(double a, const Quat &b)
{
	return b * a;
}

Quat
Quat::operator /(double rhs) const
{
	return Quat(a_ / rhs, b_ / rhs, c_ / rhs, d_ / rhs);
}

Quat
Quat::operator /(const Quat &rhs) const
{
	Quat out(*this);
	out /= rhs;
	return out;
}

Quat
operator /(double a, const Quat &b)
{
	return Quat(a, 0, 0, 0) / b;
}

bool
Quat::operator ==(const Quat &rhs) const
{
	return ((a_ == rhs.a_) && (b_ == rhs.b_) &&
	    (c_ == rhs.c_) && (d_ == rhs.d_));
}

bool
Quat::operator !=(const Quat &rhs) const
{
	return !(*this == rhs);
}

Quat
pow(const Quat &q, int n)
{
	if (n > 1) {
		int m = (n >> 1);
		Quat r = pow(q, m);
		r *= r;
		// n odd
		if (n & 1)
			r *= q;
		return r;
	}

	if (n == 1)
		return q;

	if (n == 0)
		return Quat(1, 0, 0, 0);

	// n < 0
	return pow(Quat(1, 0, 0, 0) / q, -n);
}

Quat
Quat::cross3(const Quat &v) const
{
	// Computes Euclidean cross product from the last three entries in the
	// quaternion
	return Quat(0,
	    c_ * v.d_ - d_ * v.c_,
	    d_ * v.b_ - b_ * v.d_,
	    b_ * v.c_ - c_ * v.b_);
}

Quat
cross3(const Quat &u, const Quat &v)
{
	return u.cross3(v);
}

double
Quat::dot3(const Quat &b) const
{
	// Computes Euclidean dot product from the last three entries in the
	// quaternion
	return (b_ * b.b_ + c_ * b.c_ + d_ * b.d_);
}

double
dot3(const Quat &a, const Quat &b)
{
	return a.dot3(b);
}

static G3VectorDouble
vec_abs(const G3VectorQuat &a)
{
	G3VectorDouble out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = abs(a[i]);
	return out;
}

static G3VectorQuat
vec_neg(const G3VectorQuat &a)
{
	return -1 * a;
}

static G3VectorDouble
vec_real(const G3VectorQuat &a)
{
	G3VectorDouble out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = real(a[i]);
	return out;
}

static G3Timestream
ts_abs(const G3TimestreamQuat &a)
{
	G3Timestream out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = abs(a[i]);
	return out;
}

static G3Timestream
ts_real(const G3TimestreamQuat &a)
{
	G3Timestream out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = real(a[i]);
	return out;
}

static G3TimestreamQuat
ts_neg(const G3TimestreamQuat &a)
{
	return -1 * a;
}

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
	for (auto &i: a)
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
operator /(const G3VectorQuat &a, const Quat &b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3VectorQuat
operator /(const Quat &a, const G3VectorQuat &b)
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
	for (auto &i: a)
		i /= b;
	return a;
}

G3VectorQuat &
operator /=(G3VectorQuat &a, const Quat &b)
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
operator *(const G3VectorQuat &a, const Quat &b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3VectorQuat
operator *(const Quat &b, const G3VectorQuat &a)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = b*a[i];
	return out;
}

G3VectorQuat &
operator *=(G3VectorQuat &a, const Quat &b)
{
	for (auto &i: a)
		i *= b;
	return a;
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
	for (auto &i: a)
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
operator /(const G3TimestreamQuat &a, const Quat &b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3TimestreamQuat
operator /(const Quat &a, const G3TimestreamQuat &b)
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
	for (auto &i: a)
		i /= b;
	return a;
}

G3TimestreamQuat &
operator /=(G3TimestreamQuat &a, const Quat &b)
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
operator *(const G3TimestreamQuat &a, const Quat &b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3TimestreamQuat
operator *(const Quat &b, const G3TimestreamQuat &a)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = b*a[i];
	return out;
}

G3TimestreamQuat &
operator *=(G3TimestreamQuat &a, const Quat &b)
{
	for (auto &i: a)
		i *= b;
	return a;
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


G3_SERIALIZABLE_CODE(Quat);
G3_SERIALIZABLE_CODE(G3Quat);
G3_SERIALIZABLE_CODE(G3VectorQuat);
G3_SERIALIZABLE_CODE(G3TimestreamQuat);
G3_SERIALIZABLE_CODE(G3MapQuat);
G3_SERIALIZABLE_CODE(G3MapVectorQuat);

namespace {
static int
Quat_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	py::handle<> self(py::borrowed(obj));
	py::object selfobj(self);
	py::extract<Quat &> ext(selfobj);
	if (!ext.check()) {
		PyErr_SetString(PyExc_ValueError, "Invalid quat");
		view->obj = NULL;
		return -1;
	}
	Quat &q = ext();

	view->obj = obj;
	view->buf = (void*)&q;
	view->len = 4 * sizeof(double);
	view->readonly = 0;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;

	view->ndim = 1;
	view->internal = NULL;
	view->shape = NULL;
	view->strides = NULL;
	view->suboffsets = NULL;

	Py_INCREF(obj);

	return 0;
}

static int
G3VectorQuat_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	py::handle<> self(py::borrowed(obj));
	py::object selfobj(self);
	py::extract<G3VectorQuatPtr> ext(selfobj);
	if (!ext.check()) {
		PyErr_SetString(PyExc_ValueError, "Invalid vector");
		view->obj = NULL;
		return -1;
	}
	G3VectorQuatPtr q = ext();

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

static PyBufferProcs quat_bufferprocs;
static PyBufferProcs vectorquat_bufferprocs;
static PyBufferProcs timestreamquat_bufferprocs;

static std::string
quat_str(const Quat &q)
{
	std::ostringstream oss;
	oss << q;
	return oss.str();
}

static std::string
quat_repr(const py::object &q)
{
	std::ostringstream oss;
	oss << py::extract<std::string>(q.attr("__class__").attr("__module__"))() << ".";
	oss << py::extract<std::string>(q.attr("__class__").attr("__name__"))();
	oss << py::extract<const Quat &>(q)();
	return oss.str();
}
}

std::shared_ptr<Quat>
quat_container_from_object(py::object v)
{
	// There's a chance this is actually a copy operation, so try that first
	py::extract<Quat &> extv(v);
	if (extv.check())
		return std::make_shared<Quat>(extv());

	Quat q;

	std::string format;
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_STRIDES) == -1)
		goto slowpython;

	try {
		format = check_buffer_format(view.format);
	} catch (py::buffer_error &e) {
		PyBuffer_Release(&view);
		log_fatal("%s", e.what());
	}

#define QELEM(t, i) *((t *)((char *)view.buf + i*view.strides[0]))
#define QUATI(t) Quat(QELEM(t, 0), QELEM(t, 1), QELEM(t, 2), QELEM(t, 3))

	if (view.ndim != 1 || view.shape[0] != 4) {
		PyBuffer_Release(&view);
		goto slowpython;
	} else if (format == "d") {
		q = QUATI(double);
	} else if (format == "f") {
		q = QUATI(float);
	} else if (format == "i") {
		q = QUATI(int);
	} else if (format == "l") {
		q = QUATI(long);
	} else {
		PyBuffer_Release(&view);
		goto slowpython;
	}
	PyBuffer_Release(&view);
	return std::make_shared<Quat>(q);

#undef QELEM
#undef QUATI

slowpython:
	PyErr_Clear();
	std::vector<double> xv;
	py::container_utils::extend_container(xv, v);
	if (xv.size() != 4)
		throw std::runtime_error("Invalid quat");

	return std::make_shared<Quat>(xv[0], xv[1], xv[2], xv[3]);
}

template <typename T>
std::shared_ptr<T>
quat_vec_container_from_object(py::object v)
{
	// There's a chance this is actually a copy operation, so try that first
	py::extract<T &> extv(v);
	if (extv.check())
		return std::make_shared<T>(extv());

	std::string format;
        std::shared_ptr<T> x(new (T));
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_STRIDES) == -1)
		goto slowpython;

	try {
		format = check_buffer_format(view.format);
	} catch (py::buffer_error &e) {
		PyBuffer_Release(&view);
		log_fatal("%s", e.what());
	}

#define QELEM(t, i, j) *((t *)((char *)view.buf + i*view.strides[0] + j*view.strides[1]))
#define QUATI(t, i) Quat(QELEM(t, i, 0), QELEM(t, i, 1), QELEM(t, i, 2), QELEM(t, i, 3))
#define QUATV(t) \
	for (size_t i = 0; i < (size_t) view.shape[0]; i++) \
		(*x)[i] = QUATI(t, i);

	x->resize(view.shape[0]);
	if (view.ndim != 2 || view.shape[1] != 4) {
		PyBuffer_Release(&view);
		goto slowpython;
	} else if (PyBuffer_IsContiguous(&view, 'C') &&
	    strcmp(view.format, "d") == 0 &&
	    view.strides[0] == 4*sizeof(double) &&
	    view.strides[1] == sizeof(double)) {
		// Packed and simple, use memcpy()
		memcpy((void *)&(*x)[0], view.buf, view.len);
	} else if (format == "d") {
		QUATV(double);
	} else if (format == "f") {
		QUATV(float);
	} else if (format == "i") {
		QUATV(int);
	} else if (format == "l") {
		QUATV(long);
	} else {
		PyBuffer_Release(&view);
		goto slowpython;
	}
	PyBuffer_Release(&view);
	return x;

#undef QELEM
#undef QUATI
#undef QUATV

slowpython:
	x->resize(0);
	PyErr_Clear();
	py::container_utils::extend_container(*x, v);

	return x;
}

template <>
G3VectorQuatPtr container_from_object(py::object v)
{
	return quat_vec_container_from_object<G3VectorQuat>(v);
}

template <>
G3TimestreamQuatPtr container_from_object(py::object v)
{
	return quat_vec_container_from_object<G3TimestreamQuat>(v);
}


PYBINDINGS("core", scope)
{
	auto q =
	register_class_copyable<Quat>(scope, "Quat",
	    "Representation of a quaternion. Data in a,b,c,d.")
	    .def(py::init<>())
	    .def(py::init<const Quat &>())
	    .def(py::init<double, double, double, double>(),
	        py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"),
	        "Create a quaternion from its four elements.")
	    .def("__init__", py::make_constructor(quat_container_from_object,
	        py::default_call_policies(), (py::arg("data"))),
	        "Create a quaternion from a numpy array")
	    .def(g3frameobject_picklesuite<Quat>())
	    .def_property_readonly("a", &Quat::a, "Scalar component")
	    .def_property_readonly("b", &Quat::b, "First vector component")
	    .def_property_readonly("c", &Quat::c, "Second vector component")
	    .def_property_readonly("d", &Quat::d, "Third vector component")
	    .def_property_readonly("real", &Quat::real,
	        "The real (scalar) part of the quaternion")
	    .def_property_readonly("unreal", &Quat::unreal,
	        "The unreal (vector) part of the quaternion")
	    .def(~py::self)
	    .def(-py::self)
	    .def(py::self == py::self)
	    .def(py::self != py::self)
	    .def(py::self + py::self)
	    .def(py::self += py::self)
	    .def(py::self - py::self)
	    .def(py::self -= py::self)
	    .def(py::self * py::self)
	    .def(py::self * double())
	    .def(double() * py::self)
	    .def(py::self *= py::self)
	    .def(py::self *= double())
	    .def(pow(py::self, int()))
	    .def(py::self / py::self)
	    .def(py::self / double())
	    .def(double() / py::self)
	    .def(py::self /= py::self)
	    .def(py::self /= double())
	    .def("__abs__", &Quat::abs)
	    .def("__str__", quat_str)
	    .def("__repr__", quat_repr)
	    .def("norm", &Quat::norm, "Return the Cayley norm of the quaternion")
	    .def("vnorm", &Quat::norm, "Return the Cayley norm of the "
	        "unreal (vector) part of the quaternion")
	    .def("abs", &Quat::abs, "Return the Euclidean norm of the quaternion")
	    .def("dot3", &Quat::dot3, "Dot product of last three entries")
	    .def("cross3", &Quat::cross3, "Cross product of last three entries")
	;
	PyTypeObject *qclass = (PyTypeObject *)q.ptr();
	quat_bufferprocs.bf_getbuffer = Quat_getbuffer;
	qclass->tp_as_buffer = &quat_bufferprocs;
#if PY_MAJOR_VERSION < 3
	qclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	register_frameobject<G3Quat>(scope, "G3Quat", "Serializable quaternion")
	    .def(py::init<Quat>())
	    .def_readwrite("value", &G3Quat::value)
	;

	register_vector_of<Quat>(scope, "Quat");
	auto vq =
	register_g3vector<G3VectorQuat>(scope, "G3VectorQuat",
	    "List of quaternions. Convertible to a 4xN numpy array. "
	    "Arithmetic operations on this object are fast and provide "
	    "results given proper quaternion math rather than "
	    "element-by-element numpy-ish results.")
	    .def(~py::self)
	    .def(py::self * double())
	    .def(double() * py::self)
	    .def(py::self * py::self)
	    .def(py::self * Quat())
	    .def(Quat() * py::self)
	    .def(py::self *= double())
	    .def(py::self *= Quat())
	    .def(py::self *= py::self)
	    .def(py::self / double())
	    .def(double() / py::self)
	    .def(py::self /= double())
	    .def(py::self / py::self)
	    .def(py::self /= py::self)
	    .def(py::self / Quat())
	    .def(py::self /= Quat())
	    .def(Quat() / py::self)
	    .def(pow(py::self, int()))
	    .def("__abs__", vec_abs)
	    .def("__neg__", vec_neg)
	    .def("abs", vec_abs, "Return the Euclidean norm of each quaternion")
	    .def_property_readonly("real", vec_real,
	        "Return the real (scalar) part of each quaternion")
	;
	PyTypeObject *vqclass = (PyTypeObject *)vq.ptr();
	vectorquat_bufferprocs.bf_getbuffer = G3VectorQuat_getbuffer;
	vqclass->tp_as_buffer = &vectorquat_bufferprocs;
#if PY_MAJOR_VERSION < 3
	vqclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	auto tq =
	register_frameobject<G3TimestreamQuat, G3VectorQuat>(scope, "G3TimestreamQuat",
	    "Timestream of quaternions. Identical to a G3VectorQuat except "
	    "for the addition of start and stop times.")
	    .def(py::init<>())
	    .def("__init__", py::make_constructor(container_from_object<G3TimestreamQuat>))
	    .def(py::init<const G3TimestreamQuat &>())
	    .def(~py::self)
	    .def(py::self * double())
	    .def(double() * py::self)
	    .def(py::self * G3VectorQuat())
	    .def(py::self * Quat())
	    .def(Quat() * py::self)
	    .def(py::self *= double())
	    .def(py::self *= Quat())
	    .def(py::self *= G3VectorQuat())
	    .def(py::self / double())
	    .def(double() / py::self)
	    .def(py::self /= double())
	    .def(py::self / G3VectorQuat())
	    .def(py::self /= G3VectorQuat())
	    .def(py::self / Quat())
	    .def(py::self /= Quat())
	    .def(Quat() / py::self)
	    .def(pow(py::self, int()))
	    .def("__abs__", ts_abs)
	    .def("__neg__", ts_neg)
	    .def("abs", ts_abs, "Return the Euclidean norm of each quaternion")
	    .def_property_readonly("real", ts_real,
	        "Return the real (scalar) part of each quaternion")
	    .def_readwrite("start", &G3TimestreamQuat::start,
	        "Time of the first sample in the time stream")
	    .def_readwrite("stop", &G3TimestreamQuat::stop,
	        "Time of the final sample in the timestream")
	    .def_property_readonly("sample_rate", &G3TimestreamQuat::GetSampleRate,
	        "Computed sample rate of the timestream.")
	    .def_property_readonly("n_samples", &G3TimestreamQuat::size,
	        "Number of samples in the timestream. Equivalent to len(ts)")
	;
	PyTypeObject *tqclass = (PyTypeObject *)tq.ptr();
	timestreamquat_bufferprocs.bf_getbuffer = G3VectorQuat_getbuffer;
	tqclass->tp_as_buffer = &timestreamquat_bufferprocs;
#if PY_MAJOR_VERSION < 3
	tqclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	scitbx::boost_python::container_conversions::from_python_sequence<G3TimestreamQuat, scitbx::boost_python::container_conversions::variable_capacity_policy>();
	py::implicitly_convertible<G3TimestreamQuatPtr, G3VectorQuatPtr>();
	py::implicitly_convertible<G3TimestreamQuatPtr, G3VectorQuatConstPtr>();

	register_g3map<G3MapQuat>(scope, "G3MapQuat", "Mapping from strings to "
	    "quaternions.");
	register_g3map<G3MapVectorQuat>(scope, "G3MapVectorQuat", "Mapping from "
	    "strings to lists of quaternions.");
}
