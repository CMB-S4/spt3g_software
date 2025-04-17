#ifndef _G3_PIPELINEINFO_H
#define _G3_PIPELINEINFO_H

#include <G3Frame.h>
#include <vector>
#include <map>

class G3ModuleArg : public G3FrameObject {
public:
	std::string repr;
	G3FrameObjectPtr object;

	G3ModuleArg(const std::string &r) : repr(r) {}
	G3ModuleArg(const std::string &r, G3FrameObjectPtr obj) :
	    repr(r), object(obj) {}
	G3ModuleArg() {}

	template <class A> void serialize(A &ar, unsigned v);

	std::string Description() const;

	bool operator ==(const G3ModuleArg &other) const { return other.repr == repr; };
};

class G3ModuleConfig : public G3FrameObject {
public:
	std::string modname;
	std::string instancename;

	std::map<std::string, G3ModuleArg> config;

        template <class A> void load(A &ar, unsigned v);
        template <class A> void save(A &ar, unsigned v) const;

	std::string Description() const;

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

G3_POINTERS(G3ModuleArg);
G3_POINTERS(G3ModuleConfig);
G3_POINTERS(G3PipelineInfo);

G3_SERIALIZABLE(G3ModuleArg, 1);
G3_SPLIT_SERIALIZABLE(G3ModuleConfig, 2);
G3_SERIALIZABLE(G3PipelineInfo, 2);

#endif

