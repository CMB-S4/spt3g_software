#ifndef _G3_TIMESTAMP_H
#define _G3_TIMESTAMP_H

#include <stdint.h>
#include <G3Frame.h>

typedef int64_t G3TimeStamp;

class G3Time : public G3FrameObject {
public:
	G3TimeStamp time;

	G3Time(): time(0) {}
	G3Time(int y, int d, int h, int m, int s, int ss);
	G3Time(std::string t);
	G3Time(G3TimeStamp t): time(t) {}
	G3Time(const G3Time &t): time(t.time) {}

	template <class A> void serialize(A &ar, const unsigned v);
	std::string Description() const;
	std::string isoformat() const;
	bool operator==(const G3Time & other) const;
	bool operator>(const G3Time & other) const;
	bool operator<(const G3Time & other) const;
	bool operator<=(const G3Time & other) const;
	bool operator>=(const G3Time & other) const;
	bool operator!=(const G3Time & other) const;
	G3Time &operator=(const G3Time & other) { time = other.time; return *this; }

	G3Time operator +(G3TimeStamp t) const;
	G3Time operator -(G3TimeStamp t) const;
	G3Time &operator +=(G3TimeStamp t);
	G3Time &operator -=(G3TimeStamp t);

	operator double() const;
	operator long() const;

	std::string GetFileFormatString() const;
	static G3Time Now();

	double GetMJD();
	void SetMJD(double);

private:
	SET_LOGGER("G3Time");
};

G3_POINTERS(G3Time);
G3_SERIALIZABLE(G3Time, 1);

#endif

