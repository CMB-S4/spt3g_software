#include <pybindings.h>
#include <container_pybindings.h>
#include <G3Vector.h>
#include <complex>
#include "int_storage.h"
#include <pybind11/complex.h>

G3_SPLIT_SERIALIZABLE_CODE(G3VectorInt);
G3_SERIALIZABLE_CODE(G3VectorBool);
G3_SERIALIZABLE_CODE(G3VectorDouble);
G3_SERIALIZABLE_CODE(G3VectorComplexDouble);
G3_SERIALIZABLE_CODE(G3VectorString);
G3_SERIALIZABLE_CODE(G3VectorVectorString);
G3_SERIALIZABLE_CODE(G3VectorFrameObject);
G3_SERIALIZABLE_CODE(G3VectorUnsignedChar);
G3_SERIALIZABLE_CODE(G3VectorTime);

/* Special load/save for int64_t. */

template <>
template <class A>
void G3Vector<int64_t>::load(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	int store_bits = 32;
	if (v >= 2)
		ar & cereal::make_nvp("store_bits", store_bits);

	switch(store_bits) {
	case 64:
		ar & cereal::make_nvp("vector",
		    cereal::base_class<std::vector<int64_t> >(this));
		break;
	case 32:
		load_as<A, int32_t>(ar, *this);
		break;
	case 16:
		load_as<A, int16_t>(ar, *this);
		break;
	case 8:
		load_as<A, int8_t>(ar, *this);
		break;
	}
}

template <>
template <class A>
void G3Vector<int64_t>::save(A &ar, const unsigned v) const
{
	// v == 2
	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	// Count the interesting bits, and convert to nearest power of 2.
	int sig_bits = bit_count(*this);
	int store_bits = 8;
	while (store_bits < sig_bits)
		store_bits *= 2;
	ar & cereal::make_nvp("store_bits", store_bits);
	switch(store_bits) {
	case 8:
		save_as<A, int8_t>(ar, *this);
		break;
	case 16:
		save_as<A, int16_t>(ar, *this);
		break;
	case 32:
		save_as<A, int32_t>(ar, *this);
		break;
	default:
		ar & cereal::make_nvp("vector",
		    cereal::base_class<std::vector<int64_t> >(this));
	}		
}

template <typename T, typename V>
auto fill_vector(py::buffer_info &info) {
	std::shared_ptr<V> vec(new V);
	vec->reserve((size_t) info.shape[0]);
	for (size_t i = 0; i < (size_t) info.shape[0]; i++) \
		vec->push_back(*(T *)((char *)info.ptr + i * info.strides[0]));
	return vec;
}

template <typename T, typename V>
auto fill_complex_vector(py::buffer_info &info) {
	std::shared_ptr<V> vec(new V);
	vec->resize((size_t) info.shape[0]);
	for (size_t i = 0; i < (size_t) info.shape[0]; i++) \
		(*vec)[i] = ((T *)info.ptr)[i];
	return vec;
}

template <typename V>
auto vector_from_python(const py::buffer &buf) {
	using T = typename V::value_type;

	auto info = buf.request();
	if (info.ndim != 1)
		throw py::type_error("Only valid 1D buffers can be copied to a vector");

	// handle contiguous case
	if (py::detail::compare_buffer_info<T>::compare(info)
	    && info.strides[0] == sizeof(T)
	    && info.itemsize == sizeof(T)) {
		T *p = static_cast<T *>(info.ptr);
		return std::make_shared<V>(p, p + info.shape[0]);
	}

	std::string format = check_buffer_format(info.format);

	if (format == "d")
		return fill_vector<double, V>(info);
	if (format == "f")
		return fill_vector<float, V>(info);
	if (format == "n")
		return fill_vector<ssize_t, V>(info);
	if (format == "N")
		return fill_vector<size_t, V>(info);
	if (format == "?")
		return fill_vector<bool, V>(info);
	if (format == "i")
		return fill_vector<int, V>(info);
	if (format == "I")
		return fill_vector<unsigned int, V>(info);
	if (format == "l")
		return fill_vector<long, V>(info);
	if (format == "L")
		return fill_vector<unsigned long, V>(info);
	if (format == "q")
		return fill_vector<int64_t, V>(info);
	if (format == "Q")
		return fill_vector<uint64_t, V>(info);

	throw py::type_error(std::string("Invalid buffer format ") + info.format);
}

template <typename V>
auto complex_vector_from_python(const py::buffer &buf) {
	auto info = buf.request();
	if (info.ndim != 1)
		throw py::type_error("Only valid 1D buffers can be copied to a vector");

	std::string format = check_buffer_format(info.format);

	if (format == "Zd")
		return fill_complex_vector<std::complex<double>, V>(info);
	if (format == "Zf")
		return fill_complex_vector<std::complex<float>, V>(info);

	auto vd = vector_from_python<std::vector<double> >(buf);
	return std::make_shared<V>(vd->begin(), vd->end());
}

template <typename V>
auto vector_buffer_info(V &v) {
	using T = typename V::value_type;

	return py::buffer_info(v.data(), sizeof(T),
	    py::format_descriptor<T>::format(), 1, {v.size()}, {sizeof(T)});
}

template <typename V>
auto time_vector_buffer_info(V &v) {
	using T = G3TimeStamp;

	G3Time potemkin[2];
	static ssize_t strides = (uintptr_t)&potemkin[1] - (uintptr_t)&potemkin[0];

	return py::buffer_info(&(v[0].time), sizeof(T),
	    py::format_descriptor<T>::format(), 1, {v.size()}, {strides});
}

namespace pybind11 {
	template <>
	struct format_descriptor<G3Time> : public format_descriptor<G3TimeStamp> {};
}

// specialize vector buffer implementation for time
template <typename V, typename C, typename... Args>
struct vector_buffer<G3Time, V, C, Args...> {
	static void impl(C &cls) {
		cls.def_buffer(&time_vector_buffer_info<V>);
		cls.def(py::init(&vector_from_python<V>));
		py::implicitly_convertible<py::buffer, V>();
	}
};

#define VECTOR_BUFFER_OF(T) \
template <typename V, typename C, typename... Args> \
struct vector_buffer<T, V, C, Args...> { \
	static void impl(C &cls) { \
		cls.def_buffer(&vector_buffer_info<V>); \
		cls.def(py::init(&vector_from_python<V>)); \
		py::implicitly_convertible<py::buffer, V>(); \
	} \
}

VECTOR_BUFFER_OF(float);
VECTOR_BUFFER_OF(double);
VECTOR_BUFFER_OF(int64_t);
VECTOR_BUFFER_OF(uint64_t);
VECTOR_BUFFER_OF(int32_t);
VECTOR_BUFFER_OF(uint32_t);

#define COMPLEX_VECTOR_BUFFER_OF(T) \
template <typename V, typename C, typename... Args> \
struct vector_buffer<std::complex<T>, V, C, Args...> { \
	static void impl(C &cls) { \
		cls.def_buffer(&vector_buffer_info<V>); \
		cls.def(init(&complex_vector_from_python<V>)); \
		py::implicitly_convertible<py::buffer, V>(); \
	} \
}

COMPLEX_VECTOR_BUFFER_OF(float);
COMPLEX_VECTOR_BUFFER_OF(double);


PYBINDINGS("core", scope) {

	register_vector_of<float>(scope, "Float", py::buffer_protocol());
	register_vector_of<double>(scope, "Double", py::buffer_protocol());
	register_g3vector<G3VectorDouble>(scope, "G3VectorDouble", py::buffer_protocol(),
	    "Array of floats. Treat as a serializable version of "
	    "numpy.array(dtype=float64). Can be efficiently cast to and from "
	    "numpy arrays.");

	register_vector_of<std::complex<float> >(scope, "ComplexFloat",
	    py::buffer_protocol());
	register_vector_of<std::complex<double> >(scope, "ComplexDouble",
	    py::buffer_protocol());
	register_g3vector<G3VectorComplexDouble>(scope, "G3VectorComplexDouble",
	    py::buffer_protocol(),
	    "Array of complex floats. Treat as a serializable version of "
	    "numpy.array(dtype=complex128). Can be efficiently cast to and from "
	    "numpy arrays.");

	register_vector_of<int64_t>(scope, "Int64", py::buffer_protocol());
	register_vector_of<uint64_t>(scope, "UInt64", py::buffer_protocol());
	register_vector_of<int32_t>(scope, "Int", py::buffer_protocol());
	register_vector_of<uint32_t>(scope, "UInt", py::buffer_protocol());
	register_g3vector<G3VectorInt>(scope, "G3VectorInt", py::buffer_protocol(),
	    "Array of integers. Treat as a serializable version of "
	    "numpy.array(dtype=int64). Can be efficiently cast to and from "
	    "numpy arrays.");

	register_vector_of<bool>(scope, "Bool");
	register_g3vector<G3VectorBool>(scope, "G3VectorBool", "List of booleans.");

	register_vector_of<std::string>(scope, "String");
	register_g3vector<G3VectorString>(scope, "G3VectorString", "List of strings.");

	register_vector_of<G3VectorString>(scope, "G3VectorString");
	register_g3vector<G3VectorVectorString>(scope, "G3VectorVectorString",
	    "List of lists of strings.");

	register_g3vector<G3VectorFrameObject>(scope, "G3VectorFrameObject",
	    "List of generic frame objects. Can lead to paradoxes; avoid use of "
	    "this class unless you are sure you need it.");

	register_vector_of<unsigned char>(scope, "UnsignedChar");
	register_g3vector<G3VectorUnsignedChar>(scope, "G3VectorUnsignedChar",
	    "List of 8-bit integers");

	register_vector_of<G3Time>(scope, "G3Time", py::buffer_protocol());
	register_g3vector<G3VectorTime>(scope, "G3VectorTime", py::buffer_protocol(),
	    "List of times.");
}
