#ifndef _G3_PIPELINEINFO_H
#define _G3_PIPELINEINFO_H

#include <G3Frame.h>
#include <vector>
#include <map>

class G3ModuleConfig : public G3FrameObject {
public:
	std::string modname;
	std::string instancename;

	std::map<std::string, boost::python::object> config;

        template <class A> void load(A &ar, unsigned v);
        template <class A> void save(A &ar, unsigned v) const;

	std::string Description() const;
	std::string Summary() const;

	bool operator ==(const G3ModuleConfig &) const;
};

class G3PipelineInfo : public G3FrameObject {
public:
	std::string vcs_url;
	std::string vcs_branch;
	std::string vcs_revision;
	bool vcs_localdiffs;
	std::string vcs_versionname;
	std::string vcs_githash; // If available

	std::string hostname;
	std::string user;

	std::vector<G3ModuleConfig> modules;

	template <class A> void serialize(A &ar, unsigned v);

	std::string Description() const;
	std::string Summary() const;
};

G3_POINTERS(G3PipelineInfo);

namespace cereal {
        template <class A> struct specialize<A, G3ModuleConfig, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(G3ModuleConfig, 1);
G3_SERIALIZABLE(G3PipelineInfo, 1);

#endif

