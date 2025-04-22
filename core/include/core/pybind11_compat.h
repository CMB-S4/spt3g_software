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

namespace detail {

// template magic to collect all the py::arg()'s into one object for boost
template <typename... Args>
auto extract_keywords(Args&&... args) {
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

	return std::pair(extra, std::move(kwds));
}

}

class module_ : public scope
{
public:
	static auto import(const std::string &name) {
		return py::import(name.c_str());
	}

	py::object def_submodule(const std::string &name);

	template <typename Func>
	void def(const char *name, Func && fn) {
		py::def(name, std::forward<Func>(fn));
	}

	template <typename Func>
	void def(const char *name, Func && fn, const char *doc) {
		py::def(name, std::forward<Func>(fn), doc);
	}

	template <size_t N, typename Func>
	void def(const char *name, Func && fn, py::detail::keywords<N> args, const char *doc) {
		py::def(name, std::forward<Func>(fn), args, doc);
	}

	template <typename Func, typename... Args>
	void def(const char *name, Func && fn, Args &&... args) {
		auto [extra, kwargs] = py::detail::extract_keywords(std::forward<Args>(args)...);
		std::apply([&, &kw=kwargs](auto&&... a) {
			py::def(name, std::forward<Func>(fn), kw, a...);
		}, extra);
	}
};

// class with pybind11-like functions
template <typename...T>
class compat_class_
{
public:
	typedef compat_class_<T...> self;

	compat_class_(const char *name, const char *docstring = 0)
		: cls_(name, docstring, py::no_init) {};

	auto ptr() { return cls_.ptr(); }

	template <typename... Args>
	self& add_property(Args &&... args) {
		cls_.add_property(std::forward<Args>(args)...);
		return *this;
	}

	template <typename... Args>
	self& add_static_property(Args &&... args) {
		cls_.add_static_property(std::forward<Args>(args)...);
		return *this;
	}

	template <typename... Args>
	self& def_property(const char *name, Args &&... args) {
		cls_.add_property(name, std::forward<Args>(args)...);
		return *this;
	}

	template <typename... Args>
	self& def_property_readonly(const char *name, Args &&... args) {
		cls_.add_property(name, std::forward<Args>(args)...);
		return *this;
	}

	template <typename... Args>
	self& def_readwrite(Args &&... args) {
		cls_.def_readwrite(std::forward<Args>(args)...);
		return *this;
	}

	template <typename... Args>
	self& def_readonly(Args &&... args) {
		cls_.def_readonly(std::forward<Args>(args)...);
		return *this;
	}

	self& staticmethod(const char *name) { cls_.staticmethod(name); return *this; }

	template <typename Func>
	self& def(const char *name, Func && fn) {
		cls_.def(name, std::forward<Func>(fn));
		return *this;
	}

	template <typename Func>
	self& def(const char *name, Func && fn, const char *doc) {
		cls_.def(name, std::forward<Func>(fn), doc);
		return *this;
	}

	template <typename Func, typename... Args>
	self& def(const char *name, Func && fn, Args && ... args) {
		if (sizeof...(args) == 0) {
			cls_.def(name, std::forward<Func>(fn));
			return *this;
		}
		auto [extra, kwargs] = py::detail::extract_keywords(std::forward<Args>(args)...);
		std::apply([&, &kw=kwargs](auto&&... a) {
			cls_.def(name, std::forward<Func>(fn), kw, a...);
		}, extra);
		return *this;
	}

	template <class V>
	self& def(const def_visitor<V> &visitor) {
		cls_.def(visitor);
		return *this;
	}

	template <typename... InitArgs, typename... Args>
	self& def(const py::init<InitArgs...> &initf, Args &&... args) {
		if (sizeof...(args) == 0) {
			cls_.def(initf);
			return *this;
		}
		auto [extra, kwargs] = py::detail::extract_keywords(std::forward<Args>(args)...);
		auto vis = std::apply([&, &kw=kwargs](auto&&... a) {
			return py::init<InitArgs...>(kw, a...);
		}, extra);
		cls_.def(vis);
		return *this;
	}

	template <typename... Args>
	self& def_static(const char *name, Args &&... args) {
		this->def(name, std::forward<Args>(args)...);
		cls_.staticmethod(name);
		return *this;
	}

	self& def_buffer(PyBufferProcs &p) {
		PyTypeObject *obj = (PyTypeObject *)cls_.ptr();
		obj->tp_as_buffer = &p;
		return *this;
	}

private:
      py::class_<T...> cls_;
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
