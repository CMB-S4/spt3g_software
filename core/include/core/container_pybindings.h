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
static std::shared_ptr<T>
container_from_object(py::object v)
{
	std::shared_ptr<T> x(new T());

	py::container_utils::extend_container(*x, v);
	return x;
}

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

// Register naive python conversions to and from std::vector<T>.
// See below if you also want efficient numpy conversions.
template <typename T>
py::object
register_vector_of(std::string name)
{
	name += "Vector";
	py::object cls = py::class_<std::vector<T> >(name.c_str())
	    .def("__init__", py::make_constructor(container_from_object<std::vector<T> >))
	    .def("__repr__", vec_repr<T>)
	    .def(py::vector_indexing_suite<std::vector<T>, true>())
	;
	scitbx::boost_python::container_conversions::from_python_sequence<std::vector<T>,  scitbx::boost_python::container_conversions::variable_capacity_policy>();
        return cls;
}

template <typename T, bool proxy=true>
void
register_map(std::string name, const char *docstring)
{
	py::class_<T, std::shared_ptr<T> >(name.c_str())
	    .def(py::init<const T &>())
	    .def(py::std_map_indexing_suite<T, proxy>())
	;
}

template <typename T, bool proxy=false>
void
register_g3map(std::string name, const char *docstring)
{
	// Register both regular std::map version and G3Map version,
	// with inheritance, so we can use both. Hide the std::map version
	// with _ in the Python namespace to discourage accidental use.
	register_map<std::map<typename T::key_type, typename T::mapped_type>,
	    proxy>(std::string("_") + name + "BaseMap", docstring);

	py::class_<T, py::bases<G3FrameObject,
	  std::map<typename T::key_type, typename T::mapped_type> >,
	  std::shared_ptr<T> >(name.c_str(), docstring)
	    .def(py::init<const T &>())
	    .def(py::std_map_indexing_suite<T, proxy>())
	    .def_pickle(g3frameobject_picklesuite<T>())
	;
	register_pointer_conversions<T>();
}

template <typename T>
py::class_<G3Vector<T>, py::bases<G3FrameObject, std::vector<T> >, std::shared_ptr<G3Vector<T> > >
register_g3vector(std::string name, const char *docstring)
{
	auto cls = py::class_<G3Vector<T> , py::bases<G3FrameObject, std::vector<T> >, std::shared_ptr<G3Vector<T> > >(name.c_str(), docstring) 
	    .def("__init__", py::make_constructor(container_from_object<G3Vector<T> >))
	    .def(py::vector_indexing_suite<G3Vector<T> , true>())
	    .def_pickle(g3frameobject_picklesuite<G3Vector<T> >())
	;
	scitbx::boost_python::container_conversions::from_python_sequence<G3Vector<T> ,  scitbx::boost_python::container_conversions::variable_capacity_policy>();
	register_pointer_conversions<G3Vector<T> >();

	return cls;
}

#endif
