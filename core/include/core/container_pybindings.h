#ifndef _G3_CONTAINER_PYBINDINGS_H
#define _G3_CONTAINER_PYBINDINGS_H

#include <G3Vector.h>
#include <G3Map.h>
#include <pybindings.h>
#include <serialization.h>

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/container_utils.hpp>
#include <container_conversions.h>
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include <std_map_indexing_suite.hpp>

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif

template <typename T>
std::string vec_repr(py::object self)
{
	std::stringstream s;

	s << py::extract<std::string>(self.attr("__class__").attr("__module__"))()
	    << "." << py::extract<std::string>(self.attr("__class__").attr("__name__"))()
	    << "([";

	py::extract <std::vector<T> &> extself(self);
	if (extself.check()) {
		std::vector<T> &selfobject = extself();

		int ellip_pos = -1; // Position at which to insert "..."
		if (selfobject.size() > 100)
			ellip_pos = 3;

		if (selfobject.size() > 0)
			s << selfobject[0];
		for (size_t i=1; i<selfobject.size(); ++i) {
			if ((int)i == ellip_pos) {
				s << ", ...";
				i = selfobject.size() - ellip_pos - 1;
			} else
				s << ", " << selfobject[i];
		}
	}
	s << "])";

	return s.str();
}

template <typename T>
static std::shared_ptr<T>
container_from_object(py::object v)
{
	std::shared_ptr<T> x(new T());

	py::container_utils::extend_container(*x, v);
	return x;
}


// Check buffer format for valid types
// Return format corrected for endianness
std::string check_buffer_format(std::string fmt);

// Register python conversions to and from vector types.
// Includes indexing and vector manipulation bindings
template <typename V, typename... Bases, typename... Args>
auto
register_vector(py::module_ &scope, std::string name, Args &&...args)
{
	using T = typename V::value_type;

	auto cls = register_class_copyable<V, Bases...>(scope, name, std::forward<Args>(args)...);

	cls.def("__init__", py::make_constructor(container_from_object<V>))
	    .def("__repr__", vec_repr<T>)
	    .def(py::vector_indexing_suite<V, true>())
	;
	scitbx::boost_python::container_conversions::from_python_sequence<V,
	    scitbx::boost_python::container_conversions::variable_capacity_policy>();

	return cls;
}

// Shorthand for registering the simple std::vector<T>
// Use this instead of including the pybind11/stl.h header.
template <typename T, typename... Args>
auto
register_vector_of(py::module_ &scope, std::string name, Args &&...args)
{
	name += "Vector";
	return register_vector<std::vector<T> >(scope, name, std::forward<Args>(args)...);
}

// Register a serializable vector, derived from G3FrameObject and includes
// pickling support.
template <typename V, typename... Bases, typename... Args>
auto
register_g3vector(py::module_ &scope, std::string name, Args &&...args)
{
	using U = std::vector<typename V::value_type>;

	auto cls = register_vector<V, Bases..., G3FrameObject, U>(scope,
	    name, std::forward<Args>(args)...);

	cls.def_pickle(g3frameobject_picklesuite<V>());

	register_pointer_conversions<V>();

	return cls;
}

// Register python conversions to and from map types
// Includes indexing, key/value/items views, and map manipulation bindings.
template <typename M, bool proxy=true, typename... Bases, typename... Args>
auto
register_map(py::module_ &scope, std::string name, Args &&...args)
{
	auto cls = register_class_copyable<M, Bases...>(scope, name, std::forward<Args>(args)...);

	cls.def(py::init<const M &>())
	    .def(py::std_map_indexing_suite<M, proxy>())
	;

	return cls;
}

// Shorthand for registering the single std::map<K, T>
// Use this instead of including the pybind/stl.h header.
template <typename K, typename T, bool proxy=true, typename... Args>
auto
register_map_of(py::module_ &scope, std::string name, Args &&...args)
{
	name += "Map";
	return register_map<std::map<K, T>, proxy>(scope, name, std::forward<Args>(args)...);
}

// Register a serializable map, derived from G3FrameObject and includes
// pickling support.
template <typename M, bool proxy=false, typename... Bases, typename... Args>
auto
register_g3map(py::module_ &scope, std::string name, Args &&...args)
{
	using N = std::map<typename M::key_type, typename M::mapped_type>;

	// Register both regular std::map version and G3Map version,
	// with inheritance, so we can use both. Hide the std::map version
	// with _ in the Python namespace to discourage accidental use.
	try {
		// base class may have been registered separately
		std::string base_name = std::string("_") + name + "BaseMap";
		register_map<N, proxy>(scope, base_name);
	} catch (...) {}

	auto cls = register_map<M, proxy, Bases...,  G3FrameObject, N>(scope, name,
	    std::forward<Args>(args)...);

	cls.def_pickle(g3frameobject_picklesuite<M>());

	register_pointer_conversions<M>();

	return cls;
}

#endif
