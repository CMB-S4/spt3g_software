#ifndef _G3_CONTAINER_PYBINDINGS_H
#define _G3_CONTAINER_PYBINDINGS_H

#include <pybindings.h>
#include <serialization.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif

// Bindings for vector and mapping containers

// Vector string representation
template <typename, typename, typename... Args>
void
vector_repr(const Args &...) {}

template <typename V, typename C>
auto
vector_repr(C &cl, const std::string &name)
     -> decltype(std::declval<std::ostream &>() << std::declval<typename V::value_type>(),
	void())
{
	using size_type = typename V::size_type;

	cl.def("__repr__", [name](V &v) {
		std::stringstream s;
		s << name << "([";

		int ellip_pos = -1; // Position at which to insert "..."
		if (v.size() > 100)
			ellip_pos = 3;
		if (v.size() > 0)
			s << v[0];
		for (size_type i=1; i<v.size(); ++i) {
			if ((int)i == ellip_pos) {
				s << ", ...";
				i = v.size() - ellip_pos - 1;
			} else
				s << ", " << v[i];
		}

		s << "])";
		return s.str();
	}, "Return the canonical string representation of this list.");
}

// Buffer protocol implementation
// Specialize this for particular types T (i.e. V::value_type)
// to override the default behavior.
template <typename T, typename V, typename C, typename... Args>
struct vector_buffer {
	static void impl(C &cls) {
		py::detail::vector_buffer<V, C, Args...>(cls);
	}
};

// A subclass of the pybind buffer class that ensures contiguity.
namespace PYBIND11_NAMESPACE {
class cbuffer : public buffer {
public:
	PYBIND11_OBJECT_DEFAULT(cbuffer, buffer, PyObject_CheckBuffer)

	buffer_info request_contiguous(bool writable = false) const {
		int flags = PyBUF_STRIDES | PyBUF_FORMAT | PyBUF_C_CONTIGUOUS;
		if (writable) {
			flags |= PyBUF_WRITABLE;
		}
		auto *view = new Py_buffer();
		if (PyObject_GetBuffer(m_ptr, view, flags) != 0) {
			delete view;
			throw error_already_set();
		}
		return buffer_info(view);
	}
};

namespace detail {
template <>
struct handle_type_name<cbuffer> {
	static constexpr auto name = const_name("Buffer");
};
}
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
	using C = py::class_<V, Bases..., std::shared_ptr<V> >;

	std::string rname = scope.attr("__name__").cast<std::string>() + "." + name;

	C cls(scope, name.c_str(), py::dynamic_attr(), std::forward<Args>(args)...);

	// declare buffer interface if a buffer_protocol() is passed in
	vector_buffer<typename V::value_type, V, C, Args...>::impl(cls);

	cls.def(py::init<>());

	// Register copy constructor (if possible)
	py::detail::vector_if_copy_constructible<V, C>(cls);

	// Register comparison-related operators and functions (if possible)
	py::detail::vector_if_equal_operator<V, C>(cls);

	// Register stream insertion operator (if possible)
	vector_repr<V, C>(cls, rname);

	// Modifiers require copyable vector value type
	py::detail::vector_modifiers<V, C>(cls);

	// Accessor and iterator; return by value if copyable, otherwise we return by ref + keep-alive
	py::detail::vector_accessor<V, C>(cls);

	cls.def("__bool__", [](const V &v) -> bool { return !v.empty(); },
	    "Check whether the list is nonempty");

	cls.def("__len__", [](const V &vec) { return vec.size(); });

	py::implicitly_convertible<py::iterable, V>();

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

	// Register both regular std::vector version and G3Vector version,
	// with inheritance, so we can use both. Hide the std::vector version
	// with _ in the Python namespace to discourage accidental use.
	if (!py::detail::get_type_info(typeid(U))) {
		// base class may have been registered separately
		std::string base_name = std::string("_") + name + "BaseVector";
		register_vector<U>(scope, base_name);
	}

	auto cls = register_vector<V, Bases..., U, G3FrameObject>(scope,
	    name, std::forward<Args>(args)...);

	cls.def(g3frameobject_picklesuite<V>());

	cls.def("__str__", &V::Summary)
	    .def("Summary", &V::Summary, "Short (one-line) description of the object")
	    .def("Description", &V::Description,
	        "Long-form human-readable description of the object");

	return cls;
}


// Map string representation
template <typename, typename, typename... Args>
void
map_repr(const Args &...) {}

template <typename T, typename C>
auto
map_repr(C &cls, std::string const &name)
    -> decltype(std::declval<std::ostream &>() << std::declval<typename T::key_type>()
                                               << std::declval<typename T::mapped_type>(),
       void())
{
	cls.def("__repr__", [name](T &m) {
		std::ostringstream s;
		s << name << "({";
		bool f = false;
		for (auto const &kv : m) {
			if (f)
				s << ", ";
			s << kv.first << ": " << kv.second;
			f = true;
		}
		s << "})";
		return s.str();
	}, "Return the canonical string representation of this map.");
}

// Map key/value/item iterator views for the given class
template <typename M, typename C>
void
register_map_iterators(py::module_ &scope, C &cls)
{
	using KeysView = py::detail::keys_view;
	using ValuesView = py::detail::values_view;
	using ItemsView = py::detail::items_view;

	// Wrap KeysView if it wasn't already wrapped
	if (!py::detail::get_type_info(typeid(KeysView))) {
		py::class_<KeysView> keys_view(scope, "KeysView",
		    "View of mapping keys, with length, existence and iteration operations.");
		keys_view.def("__len__", &KeysView::len);
		keys_view.def("__iter__", &KeysView::iter, py::keep_alive<0, 1>());
		keys_view.def("__contains__", &KeysView::contains);
	}
	// Similarly for ValuesView:
	if (!py::detail::get_type_info(typeid(ValuesView))) {
		py::class_<ValuesView> values_view(scope, "ValuesView",
		    "View of mapping values, with length and iteration operations.");
		values_view.def("__len__", &ValuesView::len);
		values_view.def("__iter__", &ValuesView::iter, py::keep_alive<0, 1>());
	}
	// Similarly for ItemsView:
	if (!py::detail::get_type_info(typeid(ItemsView))) {
		py::class_<ItemsView> items_view(scope, "ItemsView",
		    "View of mapping items, with length and iteration operations.");
		items_view.def("__len__", &ItemsView::len);
		items_view.def("__iter__", &ItemsView::iter, py::keep_alive<0, 1>());
	}

	cls.def("keys", [](M &m) {
		return std::unique_ptr<KeysView>(new py::detail::KeysViewImpl<M>(m));
	}, py::keep_alive<0, 1>(),
	"Return an object providing a view on the mappings's keys.");

	cls.def("values", [](M &m) {
		return std::unique_ptr<ValuesView>(new py::detail::ValuesViewImpl<M>(m));
	}, py::keep_alive<0, 1>(),
	"Return an object providing a view on the mappings's values.");

	cls.def("items", [](M &m) {
		return std::unique_ptr<ItemsView>(new py::detail::ItemsViewImpl<M>(m));
	}, py::keep_alive<0, 1>(),
	"Return an object providing a view on the mappings's items.");
}

// Register python conversions to and from map types
// Includes indexing, key/value/items views, and map manipulation bindings.
template <typename M, typename... Bases, typename... Args>
auto
register_map(py::module_ &scope, std::string name, Args &&...args)
{
	using K = typename M::key_type;
	using V = typename M::mapped_type;
	using C = py::class_<M, Bases..., std::shared_ptr<M> >;

	std::string rname = scope.attr("__name__").cast<std::string>() + "." + name;

	C cls(scope, name.c_str(), py::dynamic_attr(), std::forward<Args>(args)...);

	cls
	    .def(py::init<>())
	    .def(py::init<const M&>(), "Copy constructor")
	    .def(py::init([](const py::iterable &d) {
		auto v = std::unique_ptr<M>(new M());
		for (auto item: py::dict(d))
			v->emplace(item.first.cast<K>(), item.second.cast<V>());
		return v.release();
	    }), "Iterable constructor")
	;

	// Register stream insertion operator (if possible)
	map_repr<M, C>(cls, rname);

	// Register map iterators
	register_map_iterators<M, C>(scope, cls);

	cls.def("__iter__", [](M &m) { return py::make_key_iterator(m.begin(), m.end()); },
	    py::keep_alive<0, 1>());

	cls.def("__bool__", [](const M &m) -> bool { return !m.empty(); },
	    "Check whether the map is nonempty");

	cls.def("__getitem__", [](M &m, const K &k) -> V & {
		auto it = m.find(k);
		if (it == m.end())
			throw py::key_error();
		return it->second;
	}, py::return_value_policy::reference_internal);

	cls.def("copy", [](const M &m) {
		return std::unique_ptr<M>(new M(m)).release();
	}, "Return a shallow copy of the mapping.");

	cls.def("get", [](const M &m, const K &k, const py::object &d) {
		auto it = m.find(k);
		if (it == m.end())
			return d;
		return py::cast(it->second);
	}, py::arg("key"), py::arg("default")=py::none(),
	"Return the value for key if key is in the mapping, else default.");

	cls.def("__contains__", [](M &m, const K &k) -> bool {
		auto it = m.find(k);
		if (it == m.end())
			return false;
		return true;
	});
	// Fallback for when the object is not of the key type
	cls.def("__contains__", [](M &, const py::object &) -> bool { return false; });

	// Assignment provided only if the type is copyable
	py::detail::map_assignment<M, C>(cls);

	cls.def("update", [](py::object &m, const py::iterable &v, py::kwargs kw) {
		for (auto it: py::dict(v))
			m.attr("__setitem__")(it.first.cast<K>(), it.second.cast<V>());
		for (auto it: kw)
			m.attr("__setitem__")(it.first.cast<K>(), it.second.cast<V>());
	}, py::arg("items")=py::list(),
	"Update mapping from iterable/mapping.");

	cls.def("__delitem__", [](M &m, const K &k) {
		auto it = m.find(k);
		if (it == m.end())
			throw py::key_error();
		m.erase(it);
	});

	cls.def("pop", [](M &m, const K &k) {
		auto it = m.find(k);
		if (it == m.end())
			throw py::key_error();
		auto v = it->second;
		m.erase(it);
		return v;
	}, py::arg("key"),
	"Remove specified key and return the corresponding value. "
	"If the key is not found, raise a KeyError");

	cls.def("pop", [](M &m, const K &k, const py::object &d) {
		auto it = m.find(k);
		if (it == m.end())
			return d;
		auto v = it->second;
		m.erase(it);
		return py::cast(v);
	}, py::arg("key"), py::arg("default"),
	"Remove specified key and return the corresponding value. "
	"If the key is not found, return the default");

	cls.def("clear", [](M &m) { m.clear(); },
	    "Remove all items from the mapping.");

	// Always use a lambda in case of `using` declaration
	cls.def("__len__", [](const M &m) { return m.size(); });

	py::implicitly_convertible<py::iterable, M>();

	return cls;
}

// Shorthand for registering the single std::map<K, T>
// Use this instead of including the pybind/stl.h header.
template <typename K, typename T, typename... Args>
auto
register_map_of(py::module_ &scope, std::string name, Args &&...args)
{
	name += "Map";
	return register_map<std::map<K, T> >(scope, name, std::forward<Args>(args)...);
}

// Register a serializable map, derived from G3FrameObject and includes
// pickling support.
template <typename M, typename... Bases, typename... Args>
auto
register_g3map(py::module_ &scope, std::string name, Args &&...args)
{
	using N = std::map<typename M::key_type, typename M::mapped_type>;

	// Register both regular std::map version and G3Map version,
	// with inheritance, so we can use both. Hide the std::map version
	// with _ in the Python namespace to discourage accidental use.
	if (!py::detail::get_type_info(typeid(N))) {
		// base class may have been registered separately
		std::string base_name = std::string("_") + name + "BaseMap";
		register_map<N>(scope, base_name);
	}

	auto cls = register_map<M, Bases..., N, G3FrameObject>(scope, name,
	    std::forward<Args>(args)...);

	cls.def(g3frameobject_picklesuite<M>());

	cls.def("__str__", &M::Summary)
	    .def("Summary", &M::Summary, "Short (one-line) description of the object")
	    .def("Description", &M::Description,
	        "Long-form human-readable description of the object");

	return cls;
}

#endif
