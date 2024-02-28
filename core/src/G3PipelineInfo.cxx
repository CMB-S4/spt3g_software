#include <pybindings.h>
#include <serialization.h>
#include <G3PipelineInfo.h>
#include <std_map_indexing_suite.hpp>

#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>



template <class A> void G3ModuleConfig::save(A &ar, unsigned v) const
{
	namespace bp = boost::python;

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("modname", modname);
	ar & cereal::make_nvp("instancename", instancename);

	ar << cereal::make_nvp("size", config.size());

	G3PythonContext ctx("G3ModuleConfig::save", true);

	for (auto i : config) {
		ar << cereal::make_nvp("key", i.first);
		// Serialize frame objects (e.g. skymaps used as configs)
		// directly. Serialize random python things through repr().
		if (bp::extract<G3FrameObject>(i.second).check()) {
			G3FrameObjectConstPtr fo =
			    bp::extract<G3FrameObjectConstPtr>(i.second)();
			ar << cereal::make_nvp("frameobject", true);
			ar << cereal::make_nvp("value", fo);
		} else {
			if (!Py_IsInitialized())
				continue;

			try {
				PyObject *repr = PyObject_Repr(i.second.ptr());
				bp::handle<> reprhand(repr);
				bp::object reprobj(reprhand);
				std::string reprstr =
				    bp::extract<std::string>(reprobj);

				ar << cereal::make_nvp("frameobject", false);
				ar << cereal::make_nvp("value", reprstr);
			} catch (...) {
				log_error("Exception thrown while getting "
				    "repr() of parameter %s of module %s (%s)",
				    i.first.c_str(), instancename.c_str(),
				    modname.c_str());
				throw;
			}
		}
	}
}

template <class A> void G3ModuleConfig::load(A &ar, unsigned v)
{
	G3PythonContext ctx("G3ModuleConfig::load", true);

	namespace bp = boost::python;
	bp::object global;

	if (Py_IsInitialized()) {
		bp::object main = bp::import("__main__");
		global = main.attr("__dict__");
	}

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("modname", modname);
	ar & cereal::make_nvp("instancename", instancename);

	size_t size;
	ar >> cereal::make_nvp("size", size);

	for (size_t i = 0; i < size; i++) {
		std::string key;
		bool is_frameobject;
		ar >> cereal::make_nvp("key", key);
		ar >> cereal::make_nvp("frameobject", is_frameobject);
		
		// Frame objects (e.g. skymaps used as configs) serialized
		// directly. Random python things serialized as repr(), so
		// eval() them.
		if (is_frameobject) {
			G3FrameObjectPtr fo;
			ar >> cereal::make_nvp("value", fo);
			config[key] = boost::python::object(fo);
		} else {
			std::string repr;
			ar >> cereal::make_nvp("value", repr);

			if (!Py_IsInitialized()) {
				config[key] = boost::python::object(repr);
				continue;
			}

			bp::object obj;
			try {
//				obj = bp::eval(bp::str(repr), bp::None, bp::None);
				obj = bp::eval(repr.c_str());
			} catch (const bp::error_already_set& e) {
				obj = bp::object(repr);
				PyErr_Clear();
			}
			config[key] = obj;
		}
	}
}

std::string
G3ModuleConfig::Summary() const
{
	G3PythonContext ctx("G3ModuleConfig::Summary", true);

	std::string rv = "pipe.Add(" + modname;
	for (auto i : config) {
		std::string repr = "unknown";
		if (Py_IsInitialized())
			repr = bp::extract<std::string>(
			    i.second.attr("__repr__")());
		rv += ", " + i.first + "=" + repr;
	}

	if (instancename.size() != 0 && instancename != modname)
		rv += ", name=" + instancename;
	rv += ")";
	return rv;
}

std::string
G3ModuleConfig::Description() const
{
	return Summary();
}

bool
G3ModuleConfig::operator == (const G3ModuleConfig &b) const
{
	return (b.modname == modname) && (b.instancename == instancename) &&
	    (b.config == config);
}

template <class A> void G3PipelineInfo::serialize(A &ar, unsigned v)
{
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
	rv << " modules\n";

	return rv.str();
}

static std::string
G3PipelineInfo_repr(const G3PipelineInfo &pi)
{
	std::string rv;
	rv = "pipe = spt3g.core.G3Pipeline()\n";

	for (auto i : pi.modules) {
		rv += i.Summary();
		rv += "\n";
	}
	return rv;
}

G3ModuleConfig::~G3ModuleConfig()
{
  G3PythonContext ctx("~G3ModuleConfig",true,false);
  config.clear();
}


G3_SPLIT_SERIALIZABLE_CODE(G3ModuleConfig);
G3_SERIALIZABLE_CODE(G3PipelineInfo);

PYBINDINGS("core") {
	namespace bp = boost::python;

	register_map<std::map<std::string, boost::python::object> >(
	    "StringObjectMap", "Configuration options for a module");

	EXPORT_FRAMEOBJECT(G3ModuleConfig, init<>(), "Stored configuration of a pipeline module or segment")
	    .def_readwrite("modname", &G3ModuleConfig::modname)
	    .def_readwrite("instancename", &G3ModuleConfig::instancename)
	    .def_readwrite("config", &G3ModuleConfig::config)
	    .def("__repr__", &G3ModuleConfig::Summary)
	;
	register_pointer_conversions<G3ModuleConfig>();
	register_vector_of<G3ModuleConfig>("VectorStringObjectMap");

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
	;
	register_pointer_conversions<G3PipelineInfo>();
}

