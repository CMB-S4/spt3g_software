#pragma once

#include <boost/python.hpp>
#include <string>

#include <G3Logging.h>

namespace py = boost::python;

// Python runtime context to simplify acquiring or releasing the GIL as necessary.
// To use, simply construct the context object where necessary, e.g.
//    G3PythonContext ctx("mycontext", false);
// The context destructor will clean up after itself (releasing the GIL if acquired, and
// vice versa).  If hold_gil is true, the context will ensure the GIL is held at construction,
// and released at destruction.  If hold_gil is false, the context will save the current thread
// state and release the GIL at construction, and re-acquire it at destruction.
class G3PythonContext {
public:
	G3PythonContext(std::string name, bool hold_gil=false);
	~G3PythonContext();

private:
	std::string name_;
	bool hold_;
	PyGILState_STATE gil_;
	PyThreadState *thread_;

	SET_LOGGER("G3PythonContext");
};

// Convenience class for initializing and finalizing the Python interpreter.  This class
// will initialize the python interpeter at construction (e.g. at the beginning of a C++
// compiled program), and immediately initialize the appropriate G3PythonContext depending
// on the value of hold_gil.  At destruction, it will exit the python context and finalize
// the interpreter.  The python interpreter should be initialized only once, typically at
// the beginning of the main program.
class G3PythonInterpreter {
public:
	G3PythonInterpreter();
	~G3PythonInterpreter();

private:
	bool init_;

	SET_LOGGER("G3PythonInterpreter");
};

// pybind11 compatibility functions
namespace boost { namespace python {

// template magic to collect all the py::arg()'s into one object for boost
template <typename... Args>
auto reorder_keywords(Args&&... args) {
	using T = py::detail::keywords<1>;

	auto extra = std::tuple_cat(
		([&](auto&& arg) {
			if constexpr (!std::is_same_v<std::decay_t<decltype(arg)>, T>)
				return std::forward_as_tuple(arg);
			else
				return std::tuple<>();
		}(args))...
	);

	auto kwargs = std::tuple_cat(
		([&](auto&& arg) {
			if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, T>)
				return std::forward_as_tuple(arg);
			else
				return std::tuple<>();
		}(args))...
	);

	auto kwds = std::apply([](auto&&... kw) {
		return detail::keywords<sizeof...(kw)>{kw...};
	}, kwargs);

	return std::tuple_cat(extra, std::forward_as_tuple(kwds));
}

class module_ : public scope
{
public:
	static auto import(const std::string &name) {
		return py::import(name.c_str());
	}

	py::object def_submodule(const std::string &name);

	template <typename... Args>
	auto def(Args&&... args) {
		std::apply([](auto&&... a) { py::def(a...); },
		    reorder_keywords(args...));
	}
};

class gil_scoped_release : public G3PythonContext
{
public:
	gil_scoped_release() : G3PythonContext("release", false) {};
};

class gil_scoped_acquire : public G3PythonContext
{
public:
	gil_scoped_acquire() : G3PythonContext("acquire", true) {};
};

class scoped_interpreter : public G3PythonInterpreter {};

template <typename T>
T cast(const py::object &obj)
{
	return py::extract<T>(obj)();
}

template <typename T>
bool isinstance(const py::object &obj)
{
	return py::extract<T>(obj).check();
}

class value_error : std::exception
{
public:
	std::string msg;
	value_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class index_error : std::exception
{
public:
	std::string msg;
	index_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class type_error : std::exception
{
public:
	std::string msg;
	type_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class key_error : std::exception
{
public:
	std::string msg;
	key_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

class buffer_error : std::exception
{
public:
	std::string msg;
	buffer_error(std::string msg="") : msg(msg) {};
	const char * what() const throw() { return msg.c_str(); }
};

}}
