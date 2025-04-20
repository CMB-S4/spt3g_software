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

static std::string
quat_str(const Quat &q)
{
	std::ostringstream oss;
	oss << q;
	return oss.str();
}

static std::string
quat_repr(const Quat &q)
{
	std::ostringstream oss;
	oss << "spt3g.core.Quat" << q;
	return oss.str();
}

auto quat_buffer_info(Quat &q)
{
	return py::buffer_info(&q, sizeof(double),
	    py::format_descriptor<double>::format(), 1, {4}, {sizeof(double)});
}

template <typename T>
auto vector_quat_buffer_info(T &q)
{
	std::vector<size_t> shape{q.size(), 4};
	std::vector<size_t> strides{4 * sizeof(double), sizeof(double)};

	return py::buffer_info(&(q[0]), sizeof(double),
	    py::format_descriptor<double>::format(), 2, shape, strides);
}

namespace pybind11 {
	template <>
	struct format_descriptor<Quat> : public format_descriptor<double> {};
}

std::shared_ptr<Quat>
quat_from_python(py::buffer buf)
{
	auto info = buf.request();
	if (info.ndim != 1 || info.shape[0] != 4)
		throw py::type_error("Only valid 1D buffers can be copied to a Quat");

	std::string format = check_buffer_format(info.format);

#define QELEM(t, i) *((t *)((char *)info.ptr + i*info.strides[0]))
#define QUATI(t) std::make_shared<Quat>(QELEM(t, 0), QELEM(t, 1), QELEM(t, 2), QELEM(t, 3))

	if (format == "d")
		return QUATI(double);
	if (format == "f")
		return QUATI(float);
	if (format == "i")
		return QUATI(int);
	if (format == "l")
		return QUATI(long);

#undef QELEM
#undef QUATI

	throw py::value_error(std::string("Invalid buffer format ") + info.format);
}

std::shared_ptr<Quat>
quat_from_python_iter(py::iterable obj)
{
	if (py::len(obj) != 4)
		throw py::value_error("Invalid quat");

	auto v = obj.cast<std::vector<double> >();
	return std::make_shared<Quat>(v[0], v[1], v[2], v[3]);
}

template <typename T>
std::shared_ptr<T>
vector_quat_from_python(py::buffer &buf)
{
	auto info = buf.request();
	if (info.ndim != 2 || info.shape[1] != 4)
		throw py::type_error("Only valid 2D buffers can be copied to a Quat vector");

	std::shared_ptr<T> q(new T);
	q->resize(info.shape[0]);

	// handle contiguous case
	if (py::detail::compare_buffer_info<double>::compare(info)
	    && info.strides[0] == 4 * sizeof(double)
	    && info.strides[1] == sizeof(double)
	    && info.itemsize == sizeof(double)) {
		size_t nbytes = info.shape[0] * info.shape[1] * sizeof(double);
		memcpy((void *)&(*q)[0], info.ptr, nbytes);
		return q;
	}

	std::string format = check_buffer_format(info.format);

#define QELEM(t, i, j) *((t *)((char *)info.ptr + i*info.strides[0] + j*info.strides[1]))
#define QUATI(t, i) Quat(QELEM(t, i, 0), QELEM(t, i, 1), QELEM(t, i, 2), QELEM(t, i, 3))
#define QUATV(t) \
	for (size_t i = 0; i < (size_t) info.shape[0]; i++) \
		(*q)[i] = QUATI(t, i);

	if (format == "d") {
		QUATV(double);
	} else if (format == "f") {
		QUATV(float);
	} else if (format == "i") {
		QUATV(int);
	} else if (format == "l") {
		QUATV(long);
	} else
		throw py::value_error(std::string("Invalid buffer format :") + info.format);

#undef QELEM
#undef QUATI
#undef QUATV

	return q;
}

template <typename V, typename C, typename... Args>
struct vector_buffer<Quat, V, C, Args...> {
	static void impl(C &cls) {
		cls.def_buffer(&vector_quat_buffer_info<V>);
		cls.def(py::init(&vector_quat_from_python<V>));
		py::implicitly_convertible<py::buffer, V>();
	}
};


PYBINDINGS("core", scope)
{
	register_class<Quat>(scope, "Quat", py::buffer_protocol(),
	    "Representation of a quaternion. Data in a,b,c,d.")
	    .def(py::init<>())
	    .def(py::init<const Quat &>())
	    .def(py::init<double, double, double, double>(),
	        py::arg("a"), py::arg("b"), py::arg("c"), py::arg("d"),
	        "Create a quaternion from its four elements.")
	    .def(init(&quat_from_python), py::arg("data"),
	        "Create a quaternion from a numpy array")
	    .def(init(&quat_from_python_iter), py::arg("data"),
	        "Create a quaternion from a python iterable")
	    .def_buffer(&quat_buffer_info)
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
	    .def("__pow__", [](Quat &q, int v) { return pow(q, v); },
	        py::is_operator())
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

	register_frameobject<G3Quat>(scope, "G3Quat", "Serializable quaternion")
	    .def(py::init<Quat>())
	    .def_readwrite("value", &G3Quat::value)
	;

	register_vector_of<Quat>(scope, "Quat", py::buffer_protocol());
	register_g3vector<G3VectorQuat>(scope, "G3VectorQuat", py::buffer_protocol(),
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
	    .def("__pow__", [](G3VectorQuat &q, int v) { return pow(q, v); },
	        py::is_operator())
	    .def("__abs__", vec_abs)
	    .def("__neg__", vec_neg)
	    .def("abs", vec_abs, "Return the Euclidean norm of each quaternion")
	    .def_property_readonly("real", vec_real,
	        "Return the real (scalar) part of each quaternion")
	;

	register_g3vector<G3TimestreamQuat, G3VectorQuat>(scope, "G3TimestreamQuat",
	    py::buffer_protocol(), "Timestream of quaternions. Identical to a "
	    "G3VectorQuat except for the addition of start and stop times.")
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
	    .def("__pow__", [](G3TimestreamQuat &q, int v) { return pow(q, v); },
	        py::is_operator())
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

	py::implicitly_convertible<G3TimestreamQuat, G3VectorQuat>();

	register_g3map<G3MapQuat>(scope, "G3MapQuat", "Mapping from strings to "
	    "quaternions.");
	register_g3map<G3MapVectorQuat>(scope, "G3MapVectorQuat", "Mapping from "
	    "strings to lists of quaternions.");
}
