#ifndef _G3_PIPELINEINFO_H
#define _G3_PIPELINEINFO_H

#include <G3Frame.h>
#include <G3Map.h>
#include <vector>

class G3ModuleConfig : public G3FrameObject {
public:
	std::string modname;
	std::string instancename;

	G3MapFrameObject config;

        template <class A> void load(A &ar, unsigned v);
        template <class A> void save(A &ar, unsigned v) const;

	std::string Description() const;
	std::string Summary() const;

	bool operator ==(const G3ModuleConfig &) const;

private:
	SET_LOGGER("G3ModuleConfig");
};

class G3PipelineInfo : public G3FrameObject {
public:
	std::string vcs_url;
	std::string vcs_branch;
	std::string vcs_revision;
	bool vcs_localdiffs;
	std::string vcs_versionname;
	std::string vcs_fullversion;
	std::string vcs_githash; // If available

	std::string hostname;
	std::string user;

	std::vector<G3ModuleConfig> modules;

	template <class A> void serialize(A &ar, unsigned v);

	std::string Description() const;
	std::string Summary() const;

private:
	SET_LOGGER("G3PipelineInfo");
};

G3_POINTERS(G3PipelineInfo);

namespace cereal {
        template <class A> struct specialize<A, G3ModuleConfig, cereal::specialization::member_load_save> {};
}

G3_SERIALIZABLE(G3ModuleConfig, 2);
G3_SERIALIZABLE(G3PipelineInfo, 2);

#endif

