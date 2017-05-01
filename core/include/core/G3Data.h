#ifndef _G3_DATA_H
#define _G3_DATA_H

#include <stdint.h>
#include <G3Frame.h>

class G3Bool : public G3FrameObject {
public:
	bool value;

	G3Bool(bool val = false) : value(val) {}

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
	bool operator ==(const G3Bool & other) const {return value == other.value;}
	bool operator ==(bool other) const {return value == other;}
	bool operator !() const {return !value;}
	operator bool() const {return value;}
	bool truth() const {return value;}
};

class G3Int : public G3FrameObject {
public:
	int64_t value;

	G3Int(int64_t val = 0) : value(val) {}

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
	bool operator==(const G3Int & other) const {return value == other.value;}
};

class G3Double : public G3FrameObject {
public:
	double value;

	G3Double(double val = 0) : value(val) {}

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
	bool operator==(const G3Int & other) const {return value == other.value;}
};

class G3String : public G3FrameObject {
public:
	std::string value;
	
	G3String(const std::string &val) : value(val) {}
	G3String(const char *val = "") : value(val) {}

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
	bool operator==(const G3String & other) const {return value == other.value;}
};

G3_POINTERS(G3Bool);
G3_POINTERS(G3Int);
G3_POINTERS(G3Double);
G3_POINTERS(G3String);
G3_SERIALIZABLE(G3Bool, 1);
G3_SERIALIZABLE(G3Int, 1);
G3_SERIALIZABLE(G3Double, 1);
G3_SERIALIZABLE(G3String, 1);

#endif

