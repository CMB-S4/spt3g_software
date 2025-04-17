#pragma once

#include <G3Frame.h>
#include <G3Map.h>

#include <exception>
#include <stdint.h>

using namespace std;

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

	string Description() const;
	string Summary() const;

	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(G3TimesampleMap);
G3_SERIALIZABLE(G3TimesampleMap, 0);

class g3timesample_exception : std::exception
{
	// Exception raised when internal validity checks fail.  This will
	// also be mapped to some particular Python exception type.
public:
	std::string text;
	g3timesample_exception(std::string text) :
	    text{text} {}

	std::string msg_for_python() const throw() {
		return text;
	}
};
