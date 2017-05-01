#ifndef _G3_PYBINDINGS_H
#define _G3_PYBINDINGS_H

#include <G3.h>
#include <G3Frame.h>

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/container_utils.hpp>

#include <container_conversions.h>
#include <std_map_indexing_suite.hpp>

template <typename T>
void
register_pointer_conversions()
{
	using boost::python::implicitly_convertible;

	implicitly_convertible<boost::shared_ptr<T>, G3FrameObjectPtr>();
	implicitly_convertible<boost::shared_ptr<T>, boost::shared_ptr<const T> >();
	implicitly_convertible<boost::shared_ptr<T>, G3FrameObjectConstPtr>();
}

template <typename T>
static boost::shared_ptr<T>
container_from_object(boost::python::object v)
{
	boost::shared_ptr<T> x(new T());

	boost::python::container_utils::extend_container(*x, v);
	return x;
}

template <typename T>
std::string vec_repr(boost::python::object self)
{
	using namespace boost::python;
	std::stringstream s;

	s << extract<std::string>(self.attr("__class__").attr("__module__"))()
	    << "." << extract<std::string>(self.attr("__class__").attr("__name__"))()
	    << "([";

	std::vector<T> &selfobject = extract<std::vector<T> &>(self)();
	if (selfobject.size() == 1) {
		s << selfobject[0];
	} else if (selfobject.size() > 1){
		auto i = selfobject.begin();
		while (i != selfobject.end() - 1) {
			s << *i << ", ";
			i++;
		}
		s << *i;
	}
	s << "])";

	return s.str();
}

template <typename T>
void
register_vector_of(std::string name)
{
	name += "Vector";
	using namespace boost::python;
	class_<std::vector<T> >(name.c_str())
	    .def("__init__", make_constructor(container_from_object<std::vector<T> >))
	    .def("__repr__", vec_repr<T>)
	    .def(vector_indexing_suite<std::vector<T>, true>())
	;
	scitbx::boost_python::container_conversions::from_python_sequence<std::vector<T>,  scitbx::boost_python::container_conversions::variable_capacity_policy>();
}

template <typename T, bool proxy=true>
void
register_map(std::string name, const char *docstring)
{
	using namespace boost::python;
	class_<T, boost::shared_ptr<T> >(name.c_str())
	    .def(init<const T &>())
	    .def(std_map_indexing_suite<T, proxy>())
	;
}

class G3ModuleRegistrator {
public:
	G3ModuleRegistrator(const char *mod, void (*def)());
	static void CallRegistrarsFor(const char *mod);
};

#define EXPORT_G3MODULE(mod, T, init, docstring) \
	static void registerfunc##T() { \
		using namespace boost::python; \
		class_<T, bases<G3Module>, boost::shared_ptr<T>, \
		  boost::noncopyable>(#T, docstring, init) \
		    .def_readonly("__g3module__", true) \
		; \
	} \
	static G3ModuleRegistrator register##T(mod, registerfunc##T);

#define PYBINDINGS(mod) \
	static void ___pybindings_registerfunc(); \
	static G3ModuleRegistrator ___pybindings_register(mod, ___pybindings_registerfunc); \
	static void ___pybindings_registerfunc() 

#define EXPORT_FRAMEOBJECT(T, initf, docstring) \
	boost::python::class_<T, boost::python::bases<G3FrameObject>, boost::shared_ptr<T> >(#T, docstring, boost::python::initf) \
	    .def(boost::python::init<const T &>()) \
	    .def_pickle(g3frameobject_picklesuite<T>())

#endif

