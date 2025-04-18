#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <G3PipelineInfo.h>

namespace bp = boost::python;

template <class A> void G3ModuleArg::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("repr", repr);
	ar & cereal::make_nvp("obj", object);
}

std::string
G3ModuleArg::Description() const {
	std::string rv = "G3ModuleArg(";
	if (repr.size())
		rv += repr;
	else if (!!object)
		rv += object->Summary();
	rv += ")";
	return rv;
}

static std::string inline object_repr(bp::object obj)
{
	PyObject *repr = PyObject_Repr(obj.ptr());
	bp::handle<> reprhand(repr);
	bp::object reprobj(reprhand);
	return bp::extract<std::string>(reprobj);
}

static std::string
G3ModuleArg_repr(const G3ModuleArg &arg)
{
	// Some frame objects (e.g. G3Vectors) have repr only
	// defined properly via the python interface.
	if (!arg.repr.size() && !!arg.object)
		return object_repr(bp::object(arg.object));
	return arg.repr;
}

template <class A> void G3ModuleConfig::save(A &ar, unsigned v) const
{
	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("modname", modname);
	ar & cereal::make_nvp("instancename", instancename);
	ar & cereal::make_nvp("config", config);
}

template <class A> void G3ModuleConfig::load(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("modname", modname);
	ar & cereal::make_nvp("instancename", instancename);

	if (v > 1) {
		ar & cereal::make_nvp("config", config);
		return;
	}

	size_t size;
	ar >> cereal::make_nvp("size", size);

	for (size_t i = 0; i < size; i++) {
		std::string key;
		bool is_frameobject;
		ar >> cereal::make_nvp("key", key);
		ar >> cereal::make_nvp("frameobject", is_frameobject);
		
		// Frame objects (e.g. skymaps used as configs) serialized
		// directly. Random python things serialized as repr()
		if (is_frameobject) {
			G3FrameObjectPtr fo;
			ar >> cereal::make_nvp("value", fo);
			config[key] = G3ModuleArg("", fo);
		} else {
			std::string repr;
			ar >> cereal::make_nvp("value", repr);
			config[key] = G3ModuleArg(repr);
		}
	}
}

std::string
G3ModuleConfig_repr(const G3ModuleConfig &mc)
{
	std::string rv = "pipe.Add(" + mc.modname;
	for (auto i : mc.config)
		rv += ", " + i.first + "=" + G3ModuleArg_repr(i.second);

	if (mc.instancename.size() != 0 && mc.instancename != mc.modname)
		rv += ", name=" + mc.instancename;
	rv += ")";
	return rv;
}

std::string
G3ModuleConfig::Description() const
{
	std::ostringstream rv;
	rv << "G3ModuleConfig(" << modname;
	rv << ", " << config.size() << " arguments)";
	return rv.str();
}

bool
G3ModuleConfig::operator == (const G3ModuleConfig &b) const
{
	return (b.modname == modname) && (b.instancename == instancename) &&
	    (b.config == config);
}


static bp::object
G3ModuleConfig_get(const G3ModuleConfig &mc, std::string key)
{
	auto item = mc.config.find(key);
	if (item == mc.config.end()) {
		PyErr_SetString(PyExc_KeyError, key.c_str());
		bp::throw_error_already_set();
	}

	auto arg = item->second;
	if (!!arg.object)
		return bp::object(arg.object);

	bp::object main = bp::import("__main__");
	bp::dict global = bp::dict(main.attr("__dict__"));
	global["__main__"] = main;

	try {
		return bp::eval(bp::str(arg.repr), global, global);
	} catch (const bp::error_already_set& e) {
		PyErr_Clear();
		return bp::object(arg.repr);
	}
}

static void
G3ModuleConfig_set(G3ModuleConfig &mc, std::string key, bp::object obj)
{
	std::string repr = object_repr(obj);

	bp::extract<G3FrameObjectPtr> extobj(obj);
	if (!extobj.check()) {
		mc.config[key] = G3ModuleArg(repr);
		return;
	}

	mc.config[key] = G3ModuleArg(repr, extobj());
}

static bp::list
G3ModuleConfig_keys(const G3ModuleConfig &mc)
{
	bp::list keys;

	for (auto i: mc.config)
		keys.append(i.first);

	return keys;
}

static bp::list
G3ModuleConfig_values(const G3ModuleConfig &mc)
{
	bp::list values;

	for (auto i: mc.config)
		values.append(G3ModuleConfig_get(mc, i.first));

	return values;
}

template <class A> void G3PipelineInfo::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	using namespace cereal;

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));

	ar & make_nvp("vcs_url", vcs_url);
	ar & make_nvp("vcs_branch", vcs_branch);
	ar & make_nvp("vcs_revision", vcs_revision);
	ar & make_nvp("vcs_localdiffs", vcs_localdiffs);
	ar & make_nvp("vcs_versionname", vcs_versionname);
	ar & make_nvp("vcs_githash", vcs_githash);

	ar & make_nvp("hostname", hostname);
	ar & make_nvp("user", user);

	ar & make_nvp("modules", modules);

	if (v > 1)
		ar & make_nvp("vcs_fullversion", vcs_fullversion);
}

std::string
G3PipelineInfo::Summary() const
{
	return vcs_branch + " branch, " + ((vcs_localdiffs) ? "" : "no ") +
	    "local diffs";
}

std::string
G3PipelineInfo::Description() const
{
	std::ostringstream rv;
	rv << "Branch: " << vcs_branch <<  ", " <<
	    ((vcs_localdiffs) ? "" : "no ") << "local diffs\n";
	rv << "URL: " << vcs_url << "\n";
	rv << "Revision: " << vcs_revision << "\n";
	if (vcs_versionname.size() != 0)
		rv << "Version: " << vcs_versionname << "\n";
	if (vcs_fullversion.size() != 0)
		rv << "Full version: " << vcs_fullversion << "\n";
	rv << "Run by: " << user << " on " << hostname << "\n";

	rv << modules.size();
	rv << " modules";

	return rv.str();
}

static std::string
G3PipelineInfo_repr(const G3PipelineInfo &pi)
{
	std::string rv;
	rv = "pipe = spt3g.core.G3Pipeline()";

	for (auto i : pi.modules)
		rv += "\n" + G3ModuleConfig_repr(i);
	return rv;
}

static void
G3PipelineInfo_run(const G3PipelineInfo &pi)
{
	bp::object main = bp::import("__main__");
	bp::dict global = bp::dict(main.attr("__dict__"));
	global["__main__"] = main;

	std::string pipe = G3PipelineInfo_repr(pi);
	pipe += "\npipe.Run()";

	bp::exec(bp::str(pipe), global, global);
}

G3_SERIALIZABLE_CODE(G3ModuleArg);
G3_SPLIT_SERIALIZABLE_CODE(G3ModuleConfig);
G3_SERIALIZABLE_CODE(G3PipelineInfo);

PYBINDINGS("core") {
	EXPORT_FRAMEOBJECT(G3ModuleConfig, init<>(), "Stored configuration of a pipeline module or segment")
	    .def_readwrite("modname", &G3ModuleConfig::modname)
	    .def_readwrite("instancename", &G3ModuleConfig::instancename)
	    .def("__repr__", &G3ModuleConfig_repr)
	    .def("__getitem__", &G3ModuleConfig_get)
	    .def("__setitem__", &G3ModuleConfig_set)
	    .def("keys", &G3ModuleConfig_keys)
	    .def("values", &G3ModuleConfig_values)
	;
	register_pointer_conversions<G3ModuleConfig>();
	register_vector_of<G3ModuleConfig>("ModuleConfig");

	EXPORT_FRAMEOBJECT(G3PipelineInfo, init<>(), "Stored configuration of a pipeline, including software version information")
	    .def_readwrite("vcs_url", &G3PipelineInfo::vcs_url)
	    .def_readwrite("vcs_branch", &G3PipelineInfo::vcs_branch)
	    .def_readwrite("vcs_revision", &G3PipelineInfo::vcs_revision)
	    .def_readwrite("vcs_localdiffs", &G3PipelineInfo::vcs_localdiffs)
	    .def_readwrite("vcs_versionname", &G3PipelineInfo::vcs_versionname)
	    .def_readwrite("vcs_fullversion", &G3PipelineInfo::vcs_fullversion)
	    .def_readwrite("vcs_githash", &G3PipelineInfo::vcs_githash)
	    .def_readwrite("hostname", &G3PipelineInfo::hostname)
	    .def_readwrite("user", &G3PipelineInfo::user)
	    .def_readwrite("modules", &G3PipelineInfo::modules)
	    .def("__repr__", &G3PipelineInfo_repr)
	    .def("Run", &G3PipelineInfo_run)
	;
	register_pointer_conversions<G3PipelineInfo>();
}

