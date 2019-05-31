#ifndef _G3_CONTAINER_PYBINDINGS_H
#define _G3_CONTAINER_PYBINDINGS_H

#include <G3Vector.h>
#include <G3Map.h>
#include <pybindings.h>
#include <serialization.h>
#include <std_map_indexing_suite.hpp>

template <typename T, bool proxy=false>
void
register_g3map(std::string name, const char *docstring)
{
	using namespace boost::python;

	// Register both regular std::map version and G3Map version,
	// with inheritance, so we can use both. Hide the std::map version
	// with _ in the Python namespace to discourage accidental use.
	register_map<std::map<typename T::key_type, typename T::mapped_type>,
	    proxy>(std::string("_") + name + "BaseMap", docstring);

	class_<T, bases<G3FrameObject,
	  std::map<typename T::key_type, typename T::mapped_type> >,
	  boost::shared_ptr<T> >(name.c_str(), docstring)
	    .def(init<const T &>())
	    .def(std_map_indexing_suite<T, proxy>())
	    .def_pickle(g3frameobject_picklesuite<T>())
	;
	register_pointer_conversions<T>();
}

template <typename T>
boost::python::class_<G3Vector<T>, boost::python::bases<G3FrameObject, std::vector<T> >, boost::shared_ptr<G3Vector<T> > >
register_g3vector(std::string name, const char *docstring)
{
	using namespace boost::python;
	auto cls = class_<G3Vector<T> , bases<G3FrameObject, std::vector<T> >, boost::shared_ptr<G3Vector<T> > >(name.c_str(), docstring) 
	    .def("__init__", make_constructor(container_from_object<G3Vector<T> >))
	    .def(vector_indexing_suite<G3Vector<T> , true>())
	    .def_pickle(g3frameobject_picklesuite<G3Vector<T> >())
	;
	scitbx::boost_python::container_conversions::from_python_sequence<G3Vector<T> ,  scitbx::boost_python::container_conversions::variable_capacity_policy>();
	register_pointer_conversions<G3Vector<T> >();

	return cls;
}

#endif
