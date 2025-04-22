#pragma once

#include <G3Frame.h>
#include <G3Map.h>

class G3TimesampleMap : public G3MapFrameObject {
	// Storage for a set of co-sampled data vectors, and the single
	// vector of associated timestamps.  The vectors can have
	// different types, but must be one of the explicitly handled
	// types.  Run Check() to confirm your object has been populated
	// properly.
public:
	G3VectorTime times;
	G3TimesampleMap Concatenate(const G3TimesampleMap &other) const;
	bool Check() const;
	void Sort();

	std::string Description() const;
	std::string Summary() const;

	template <class A> void serialize(A &ar, unsigned v);

private:
	SET_LOGGER("G3TimesampleMap");
};

G3_POINTERS(G3TimesampleMap);
G3_SERIALIZABLE(G3TimesampleMap, 0);
