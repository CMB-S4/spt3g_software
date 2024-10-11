#include <pybindings.h>
#include <container_pybindings.h>
#include <G3Quat.h>
#include <G3Map.h>
#include <G3Units.h>

// Quaternion utilities

std::string
G3Quat::Description() const
{
	std::ostringstream desc;
	desc << "[" << buf_[0] << ", " << buf_[1] << ", " << buf_[2] << ", " << buf_[3] << "]";
	if (versor_)
		desc << ", versor=True";
	return desc.str();
}

template <class A>
void G3Quat::serialize(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("a", buf_[0]);
	ar & cereal::make_nvp("b", buf_[1]);
	ar & cereal::make_nvp("c", buf_[2]);
	ar & cereal::make_nvp("d", buf_[3]);
	ar & cereal::make_nvp("versor", versor_);
}

typedef struct {
	double a, b, c, d;
} V1Quat;

template<class A>
void serialize(A &ar, V1Quat &q, unsigned version)
{
	using namespace cereal;
	ar & make_nvp("a", q.a);
	ar & make_nvp("b", q.b);
	ar & make_nvp("c", q.c);
	ar & make_nvp("d", q.d);
}

template<>
template <class A>
void G3VectorQuat::load(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));

	if (v > 1) {
		ar & cereal::make_nvp("vector",
		    cereal::base_class<std::vector<G3Quat> >(this));
	} else {
		std::vector<V1Quat> vec;
		ar & cereal::make_nvp("vector", vec);

		this->resize(0);
		for (auto &v: vec)
			this->push_back(G3Quat(v.a, v.b, v.c, v.d));
	}
}

template<>
template <class A>
void G3VectorQuat::save(A &ar, unsigned v) const
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("vector",
	    cereal::base_class<std::vector<G3Quat> >(this));
}

template<>
template <class A>
void G3MapQuat::load(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));

	if (v > 1) {
		ar & cereal::make_nvp("map",
		    cereal::base_class<std::map<std::string, G3Quat> >(this));
	} else {
		std::map<std::string, V1Quat> m;
		ar & cereal::make_nvp("map", m);

		this->clear();
		for (auto &i: m)
			(*this)[i.first] = G3Quat(i.second.a, i.second.b,
			    i.second.c, i.second.d);
	}
}

template<>
template <class A>
void G3MapQuat::save(A &ar, unsigned v) const
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("map",
	    cereal::base_class<std::map<std::string, G3Quat> >(this));
}

void G3Quat::versor_inplace()
{
	if (!versor_) {
		double n = norm();
		if (fabs(n - 1.0) > 1e-6)
			*this /= sqrt(n);
		versor_ = true;
	}
}

G3Quat
G3Quat::versor() const
{
	G3Quat out(*this);
	out.versor_inplace();
	return out;
}

double
G3Quat::real() const
{
	return buf_[0];
}

G3Quat
G3Quat::unreal() const
{
	if (!buf_[0])
		return *this;
	return G3Quat(0, buf_[1], buf_[2], buf_[3]);
}

G3Quat
G3Quat::conj() const
{
	return G3Quat(buf_[0], -buf_[1], -buf_[2], -buf_[3], versor_);
}

double
G3Quat::norm() const
{
	return buf_[0] * buf_[0] + buf_[1] * buf_[1] +
	    buf_[2] * buf_[2] + buf_[3] * buf_[3];
}

double
G3Quat::abs() const
{
	return sqrt(norm());
}

void *
G3Quat::buffer()
{
	return (void *)(&(buf_[0]));
}

G3Quat
G3Quat::operator ~() const
{
	return conj();
}

G3Quat &
G3Quat::operator +=(const G3Quat &rhs)
{
	buf_[0] += rhs.buf_[0];
	buf_[1] += rhs.buf_[1];
	buf_[2] += rhs.buf_[2];
	buf_[3] += rhs.buf_[3];
	versor_ = false;
	return *this;
}

G3Quat &
G3Quat::operator -=(const G3Quat &rhs)
{
	buf_[0] -= rhs.buf_[0];
	buf_[1] -= rhs.buf_[1];
	buf_[2] -= rhs.buf_[2];
	buf_[3] -= rhs.buf_[3];
	versor_ = false;
	return *this;
}

G3Quat &
G3Quat::operator *=(double rhs)
{
	buf_[0] *= rhs;
	buf_[1] *= rhs;
	buf_[2] *= rhs;
	buf_[3] *= rhs;
	versor_ = false;
	return *this;
}

G3Quat &
G3Quat::operator *=(const G3Quat &rhs)
{
	const double *vr = (const double *)(&(rhs.buf_[0]));
	double a = buf_[0] * vr[0] - buf_[1] * vr[1] - buf_[2] * vr[2] - buf_[3] * vr[3];
	double b = buf_[0] * vr[1] + buf_[1] * vr[0] + buf_[2] * vr[3] - buf_[3] * vr[2];
	double c = buf_[0] * vr[2] - buf_[1] * vr[3] + buf_[2] * vr[0] + buf_[3] * vr[1];
	double d = buf_[0] * vr[3] + buf_[1] * vr[2] - buf_[2] * vr[1] + buf_[3] * vr[0];
	buf_[0] = a;
	buf_[1] = b;
	buf_[2] = c;
	buf_[3] = d;
	if (is_versor() && rhs.is_versor())
		versor_inplace();
	else
		versor_ = false;
	return *this;
}

G3Quat &
G3Quat::operator /=(double rhs)
{
	buf_[0] /= rhs;
	buf_[1] /= rhs;
	buf_[2] /= rhs;
	buf_[3] /= rhs;
	versor_ = false;
	return *this;
}

G3Quat &
G3Quat::operator /=(const G3Quat &rhs)
{
	double n = rhs.norm();
	const double *vr = (const double *)(&(rhs.buf_[0]));
	double a =  buf_[0] * vr[0] + buf_[1] * vr[1] + buf_[2] * vr[2] + buf_[3] * vr[3];
	double b = -buf_[0] * vr[1] + buf_[1] * vr[0] - buf_[2] * vr[3] + buf_[3] * vr[2];
	double c = -buf_[0] * vr[2] + buf_[1] * vr[3] + buf_[2] * vr[0] - buf_[3] * vr[1];
	double d = -buf_[0] * vr[3] - buf_[1] * vr[2] + buf_[2] * vr[1] + buf_[3] * vr[0];
	buf_[0] = a / n;
	buf_[1] = b / n;
	buf_[2] = c / n;
	buf_[3] = d / n;
	if (is_versor() && rhs.is_versor())
		versor_inplace();
	else
		versor_ = false;
	return *this;
}

G3Quat
G3Quat::operator +(const G3Quat &rhs) const
{
	return G3Quat(buf_[0] + rhs.buf_[0], buf_[1] + rhs.buf_[1],
	    buf_[2] + rhs.buf_[2], buf_[3] + rhs.buf_[3]);
}

G3Quat
G3Quat::operator -(const G3Quat &rhs) const
{
	return G3Quat(buf_[0] - rhs.buf_[0], buf_[1] - rhs.buf_[1],
	    buf_[2] - rhs.buf_[2], buf_[3] - rhs.buf_[3]);
}

G3Quat
G3Quat::operator *(double rhs) const
{
	return G3Quat(buf_[0] * rhs, buf_[1] * rhs, buf_[2] * rhs, buf_[3] * rhs);
}

G3Quat
G3Quat::operator *(const G3Quat &rhs) const
{
	G3Quat out(*this);
	out *= rhs;
	return out;
}

G3Quat
operator *(double a, const G3Quat &b)
{
	return b * a;
}

G3Quat
G3Quat::operator /(double rhs) const
{
	return G3Quat(buf_[0] / rhs, buf_[1] / rhs, buf_[2] / rhs, buf_[3] / rhs);
}

G3Quat
G3Quat::operator /(const G3Quat &rhs) const
{
	G3Quat out(*this);
	out /= rhs;
	return out;
}

G3Quat
operator /(double a, const G3Quat &b)
{
	return G3Quat(a, 0, 0, 0) / b;
}

bool
G3Quat::operator ==(const G3Quat &rhs) const
{
	return ((buf_[0] == rhs.buf_[0]) && (buf_[1] == rhs.buf_[1]) &&
	    (buf_[2] == rhs.buf_[2]) && (buf_[3] == rhs.buf_[3]));
}

bool
G3Quat::operator !=(const G3Quat &rhs) const
{
	return !(*this == rhs);
}

G3Quat
pow(const G3Quat &q, int n)
{
	if (n > 1) {
		int m = (n >> 1);
		G3Quat r = pow(q, m);
		r *= r;
		// n odd
		if (n & 1)
			r *= q;
		return r;
	}

	if (n == 1)
		return q;

	if (n == 0)
		return G3Quat(1, 0, 0, 0);

	// n < 0
	return pow(G3Quat(1, 0, 0, 0) / q, -n);
}

G3Quat
cross3(const G3Quat &u, const G3Quat &v)
{
	// Computes Euclidean cross product from the last three entries in the
	// quaternion
	G3Quat out(0,
	    u.c()*v.d() - (u.d()*v.c()),
	    u.d()*v.b() - (u.b()*v.d()),
	    u.b()*v.c() - (u.c()*v.b()));
	if (u.is_versor() && v.is_versor())
		return out.versor();
	return out;
}

double
dot3(const G3Quat &a, const G3Quat &b)
{
	// Computes Euclidean dot product from the last three entries in the
	// quaternion
	return (a.b()*b.b() +
		a.c()*b.c() +
		a.d()*b.d());
}

static G3VectorDouble
vec_abs(const G3VectorQuat &a)
{
	G3VectorDouble out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = abs(a[i]);
	return out;
}

static G3VectorDouble
vec_real(const G3VectorQuat &a)
{
	G3VectorDouble out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = real(a[i]);
	return out;
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
	for (G3Quat &i: a)
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
operator /(const G3VectorQuat &a, const G3Quat &b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3VectorQuat
operator /(const G3Quat &a, const G3VectorQuat &b)
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
	for (G3Quat &i: a)
		i /= b;
	return a;
}

G3VectorQuat &
operator /=(G3VectorQuat &a, const G3Quat &b)
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
operator *(const G3VectorQuat &a, const G3Quat &b)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3VectorQuat
operator *(const G3Quat &b, const G3VectorQuat &a)
{
	G3VectorQuat out(a.size());
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = b*a[i];
	return out;
}

G3VectorQuat &
operator *=(G3VectorQuat &a, const G3Quat &b)
{
	for (G3Quat &i: a)
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
	for (G3Quat &i: a)
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
operator /(const G3TimestreamQuat &a, const G3Quat &b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]/b;
	return out;
}

G3TimestreamQuat
operator /(const G3Quat &a, const G3TimestreamQuat &b)
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
	for (G3Quat &i: a)
		i /= b;
	return a;
}

G3TimestreamQuat &
operator /=(G3TimestreamQuat &a, const G3Quat &b)
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
operator *(const G3TimestreamQuat &a, const G3Quat &b)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = a[i]*b;
	return out;
}

G3TimestreamQuat
operator *(const G3Quat &b, const G3TimestreamQuat &a)
{
	G3TimestreamQuat out(a.size());
	out.start = a.start; out.stop = a.stop;
	for (unsigned i = 0; i < a.size(); i++)
		out[i] = b*a[i];
	return out;
}

G3TimestreamQuat &
operator *=(G3TimestreamQuat &a, const G3Quat &b)
{
	for (G3Quat &i: a)
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


G3_SERIALIZABLE_CODE(G3Quat);
G3_SPLIT_SERIALIZABLE_CODE(G3VectorQuat);
G3_SERIALIZABLE_CODE(G3TimestreamQuat);
G3_SPLIT_SERIALIZABLE_CODE(G3MapQuat);
G3_SERIALIZABLE_CODE(G3MapVectorQuat);

namespace {
static int
G3Quat_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	bp::handle<> self(bp::borrowed(obj));
	bp::object selfobj(self);
	bp::extract<G3QuatPtr> ext(selfobj);
	if (!ext.check()) {
		PyErr_SetString(PyExc_ValueError, "Invalid quat");
		view->obj = NULL;
		return -1;
	}
	G3QuatPtr q = ext();

	view->obj = obj;
	view->buf = q->buffer();
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

	bp::handle<> self(bp::borrowed(obj));
	bp::object selfobj(self);
	bp::extract<G3VectorQuatPtr> ext(selfobj);
	if (!ext.check()) {
		PyErr_SetString(PyExc_ValueError, "Invalid vector");
		view->obj = NULL;
		return -1;
	}
	G3VectorQuatPtr q = ext();

	G3Quat potemkin[2];
	static Py_ssize_t stride0 = (uintptr_t)potemkin[1].buffer() -
	    (uintptr_t)potemkin[0].buffer();

	view->obj = obj;
	view->buf = (*q)[0].buffer();
	view->len = q->size() * stride0;
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
	view->strides[0] = stride0;
	view->strides[1] = view->itemsize;

	view->suboffsets = NULL;

	Py_INCREF(obj);

	return 0;
}

static PyBufferProcs quat_bufferprocs;
static PyBufferProcs vectorquat_bufferprocs;
static PyBufferProcs timestreamquat_bufferprocs;

static std::string
quat_repr(const G3Quat &q)
{
	std::ostringstream oss;
	oss << "spt3g.core.G3Quat(" << q.Description() << ")";
	return oss.str();
}
}

boost::shared_ptr<G3Quat>
quat_container_from_object(boost::python::object v, bool versor)
{
	// There's a chance this is actually a copy operation, so try that first
	bp::extract<G3Quat &> extv(v);
	if (extv.check())
		return boost::make_shared<G3Quat>(extv());

	boost::shared_ptr<G3Quat> x(new G3Quat(0, 0, 0, 0, versor));
	double *data = (double *)(x->buffer());

	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_STRIDES) == -1)
		goto slowpython;

#define QELEM(t, i) *((t *)((char *)view.buf + i*view.strides[0]))

	if (view.ndim != 1 || view.shape[0] != 4) {
		PyBuffer_Release(&view);
		goto slowpython;
	} else if (PyBuffer_IsContiguous(&view, 'C') &&
	    strcmp(view.format, "d") == 0 &&
	    view.strides[0] == sizeof(double)) {
		// Packed and simple, use memcpy()
		memcpy(data, (double *)((char *)view.buf), 4 * sizeof(double));
	} else if (strcmp(view.format, "d") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			data[i] = QELEM(double, i);
	} else if (strcmp(view.format, "f") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			data[i] = QELEM(float, i);
	} else if (strcmp(view.format, "i") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			data[i] = QELEM(int, i);
	} else if (strcmp(view.format, "l") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			data[i] = QELEM(long, i);
	} else {
		PyBuffer_Release(&view);
		goto slowpython;
	}
	PyBuffer_Release(&view);
	return x;

#undef QELEM

slowpython:
	PyErr_Clear();
	std::vector<double> xv;
	boost::python::container_utils::extend_container(xv, v);
	if (xv.size() != 4)
		throw std::runtime_error("Invalid quat");

	memcpy(data, &(xv[0]), 4 * sizeof(double));
	return x;
}

template <typename T>
boost::shared_ptr<T>
quat_vec_container_from_object(boost::python::object v)
{
	// There's a chance this is actually a copy operation, so try that first
	bp::extract<T &> extv(v);
	if (extv.check())
		return boost::make_shared<T>(extv());

        boost::shared_ptr<T> x(new (T));
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_STRIDES) == -1)
		goto slowpython;

#define QELEM(t, i, j) *((t *)((char *)view.buf + i*view.strides[0] + j*view.strides[1]))
#define QUATI(t, i) G3Quat(QELEM(t, i, 0), QELEM(t, i, 1), QELEM(t, i, 2), QELEM(t, i, 3))

	x->resize(view.shape[0]);
	if (view.ndim != 2 || view.shape[1] != 4) {
		PyBuffer_Release(&view);
		goto slowpython;
	} else if (PyBuffer_IsContiguous(&view, 'C') &&
	    strcmp(view.format, "d") == 0 &&
	    view.strides[0] == 4*sizeof(double) &&
	    view.strides[1] == sizeof(double)) {
		// Packed and simple, use memcpy()
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			memcpy((*x)[i].buffer(),
			    (double *)((char *)view.buf + i*view.strides[0]),
			    4 * sizeof(double));
	} else if (strcmp(view.format, "d") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			(*x)[i] = QUATI(double, i);
	} else if (strcmp(view.format, "f") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			(*x)[i] = QUATI(float, i);
	} else if (strcmp(view.format, "i") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			(*x)[i] = QUATI(int, i);
	} else if (strcmp(view.format, "l") == 0) {
		for (size_t i = 0; i < (size_t)view.shape[0]; i++)
			(*x)[i] = QUATI(long, i);
	} else {
		PyBuffer_Release(&view);
		goto slowpython;
	}
	PyBuffer_Release(&view);
	return x;

#undef QUATI
#undef QELEM

slowpython:
	x->resize(0);
	PyErr_Clear();
	boost::python::container_utils::extend_container(*x, v);

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

	object q = EXPORT_FRAMEOBJECT(G3Quat, init<>(), "Representation of a quaternion. Data in a,b,c,d.")
	     .def(init<double, double, double, double, bool>(
	         "Create a quaternion (or a unit quaternion if versor is True) from its four elements.",
	         (arg("a"), arg("b"), arg("c"), arg("d"), arg("versor")=false)))
	     .def("__init__", make_constructor(quat_container_from_object, default_call_policies(),
	         (arg("data"), arg("versor")=false)), "Create a quaternion (or versor) from a numpy array")
	     .add_property("a", &G3Quat::a, "Scalar component")
	     .add_property("b", &G3Quat::b, "First vector component")
	     .add_property("c", &G3Quat::c, "Second vector component")
	     .add_property("d", &G3Quat::d, "Third vector component")
	     .add_property("is_versor", &G3Quat::is_versor, "True if this is a unit quaternion")
	     .add_property("real", &G3Quat::real, "The real (scalar) part of the quaternion")
	     .add_property("unreal", &G3Quat::unreal, "The unreal (vector) part of the quaternion")
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
	     .def(pow(self, int()))
	     .def(self / self)
	     .def(self / double())
	     .def(double() / self)
	     .def(self /= self)
	     .def(self /= double())
	     .def("__abs__", &G3Quat::abs)
	     .def("__repr__", quat_repr)
	     .def("norm", &G3Quat::norm, "Return the Cayley norm of the quaternion")
	     .def("abs", &G3Quat::abs, "Return the Euclidean norm of the quaternion")
	     .def("versor", &G3Quat::versor, "Return a versor (unit quaternion) with the same orientation")
	     .def("dot3", dot3, "Dot product of last three entries")
	     .def("cross3", cross3, "Cross product of last three entries")
	;
	register_pointer_conversions<G3Quat>();
	PyTypeObject *qclass = (PyTypeObject *)q.ptr();
	quat_bufferprocs.bf_getbuffer = G3Quat_getbuffer;
	qclass->tp_as_buffer = &quat_bufferprocs;
#if PY_MAJOR_VERSION < 3
	qclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	register_vector_of<G3Quat>("Quat");
	object vq =
	    register_g3vector<G3Quat>("G3VectorQuat",
	     "List of quaternions. Convertible to a 4xN numpy array. "
	     "Arithmetic operations on this object are fast and provide "
	     "results given proper quaternion math rather than "
	     "element-by-element numpy-ish results.")
	     .def(~self)
	     .def(self * double())
	     .def(double() * self)
	     .def(self * self)
	     .def(self * G3Quat())
	     .def(G3Quat() * self)
	     .def(self *= double())
	     .def(self *= G3Quat())
	     .def(self *= self)
	     .def(self / double())
	     .def(double() / self)
	     .def(self /= double())
	     .def(self / self)
	     .def(self /= self)
	     .def(self / G3Quat())
	     .def(self /= G3Quat())
	     .def(G3Quat() / self)
	     .def(pow(self, int()))
	     .def("__abs__", vec_abs)
	     .def("abs", vec_abs, "Return the Euclidean norm of each quaternion")
	     .add_property("real", vec_real, "Return the real (scalar) part of each quaternion");
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
	     .def(self * G3Quat())
	     .def(G3Quat() * self)
	     .def(self *= double())
	     .def(self *= G3Quat())
	     .def(self *= G3VectorQuat())
	     .def(self / double())
	     .def(double() / self)
	     .def(self /= double())
	     .def(self / G3VectorQuat())
	     .def(self /= G3VectorQuat())
	     .def(self / G3Quat())
	     .def(self /= G3Quat())
	     .def(G3Quat() / self)
	     .def(pow(self, int()))
	     .def("__abs__", vec_abs)
	     .def("abs", vec_abs, "Return the Euclidean norm of each quaternion")
	     .add_property("real", vec_real, "Return the real (scalar) part of each quaternion")
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

	register_g3map<G3MapQuat>("G3MapQuat", "Mapping from strings to "
	    "quaternions.");
	register_g3map<G3MapVectorQuat>("G3MapVectorQuat", "Mapping from "
	    "strings to lists of quaternions.");
}
