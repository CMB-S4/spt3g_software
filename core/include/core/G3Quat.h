#ifndef _CORE_G3QUAT_H
#define _CORE_G3QUAT_H

#include <G3Vector.h>
#include <G3Map.h>

class Quat
{
public:
	Quat() : a_(0), b_(0), c_(0), d_(0) {}
	Quat(double a, double b, double c, double d) :
	    a_(a), b_(b), c_(c), d_(d) {}
	Quat(const Quat &q) : a_(q.a_), b_(q.b_), c_(q.c_), d_(q.d_) {}

	double a() const { return a_; }
	double b() const { return b_; }
	double c() const { return c_; }
	double d() const { return d_; }

	double real() const;
	Quat unreal() const;
	Quat conj() const;
	double norm() const;
	double vnorm() const;
	double abs() const;
	double dot3(const Quat &b) const;
	Quat cross3(const Quat &b) const;

	Quat operator -() const;
	Quat operator ~() const;

	Quat &operator +=(const Quat &);
	Quat &operator -=(const Quat &);
	Quat &operator *=(double);
	Quat &operator *=(const Quat &);
	Quat &operator /=(double);
	Quat &operator /=(const Quat &);

	Quat operator +(const Quat &) const;
	Quat operator -(const Quat &) const;
	Quat operator *(double) const;
	Quat operator *(const Quat &) const;
	Quat operator /(double) const;
	Quat operator /(const Quat &) const;

	bool operator ==(const Quat &) const;
	bool operator !=(const Quat &) const;

	template <class A> void serialize(A &ar, unsigned v);
private:
	double a_, b_, c_, d_;
};

std::ostream& operator<<(std::ostream& os, const Quat &);

namespace cereal {
	template <class A> struct specialize<A, Quat, cereal::specialization::member_serialize> {};
}

CEREAL_CLASS_VERSION(Quat, 1);

Quat operator *(double, const Quat &);
Quat operator /(double, const Quat &);

inline double real(const Quat &q) { return q.real(); };
inline Quat unreal(const Quat &q) { return q.unreal(); };
inline Quat conj(const Quat &q) { return q.conj(); };
inline double norm(const Quat &q) { return q.norm(); }
inline double vnorm(const Quat &q) { return q.vnorm(); }
inline double abs(const Quat &q) { return q.abs(); }

Quat pow(const Quat &, int);

Quat cross3(const Quat &a, const Quat &b);
double dot3(const Quat &a, const Quat &b);

// Frame object data wrapper

class G3Quat : public G3FrameObject {
public:
	Quat value;

	G3Quat() {}
	G3Quat(const Quat &val) : value(val) {}

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
	bool operator==(const G3Quat & other) const {return value == other.value;}
};

G3_POINTERS(G3Quat);
G3_SERIALIZABLE(G3Quat, 1);

G3VECTOR_OF(Quat, G3VectorQuat);

class G3TimestreamQuat : public G3VectorQuat
{
public:
	G3TimestreamQuat() : G3VectorQuat() {}
        G3TimestreamQuat(std::vector<Quat>::size_type s) : G3VectorQuat(s) {}
        G3TimestreamQuat(std::vector<Quat>::size_type s,
            const Quat &val) : G3VectorQuat(s, val) {}
        G3TimestreamQuat(const G3TimestreamQuat &r) : G3VectorQuat(r),
            start(r.start), stop(r.stop) {}
        G3TimestreamQuat(const G3VectorQuat &r) : G3VectorQuat(r) {}
        template <typename Iterator> G3TimestreamQuat(Iterator l, Iterator r) :
            G3VectorQuat(l, r) {}

	G3Time start, stop;
	double GetSampleRate() const;

	template <class A> void serialize(A &ar, unsigned v);

	std::string Description() const;
	std::string Summary() const { return Description(); };
};

namespace cereal {
	template <class A> struct specialize<A, G3TimestreamQuat, cereal::specialization::member_serialize> {};
}

G3_POINTERS(G3TimestreamQuat);
G3_SERIALIZABLE(G3TimestreamQuat, 1);

G3VectorQuat operator ~ (const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, double);
G3VectorQuat &operator *= (G3VectorQuat &, double);
G3VectorQuat operator / (const G3VectorQuat &, double);
G3VectorQuat operator / (double, const G3VectorQuat &);
G3VectorQuat operator / (const G3VectorQuat &, const Quat &);
G3VectorQuat operator / (const Quat &, const G3VectorQuat &);
G3VectorQuat operator / (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator /= (G3VectorQuat &, double);
G3VectorQuat &operator /= (G3VectorQuat &, const Quat &);
G3VectorQuat &operator /= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (double, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, const Quat &);
G3VectorQuat operator * (const Quat &, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, const Quat &);

G3VectorQuat pow(const G3VectorQuat &a, int b);

G3TimestreamQuat operator ~ (const G3TimestreamQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, double);
G3TimestreamQuat operator * (double, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, double);
G3TimestreamQuat operator / (double, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, const Quat &);
G3TimestreamQuat operator / (const Quat &, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, double);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, const Quat &);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat &operator *= (G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat operator * (double, const G3TimestreamQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, const Quat &);
G3TimestreamQuat operator * (const Quat &, const G3TimestreamQuat &);
G3TimestreamQuat &operator *= (G3TimestreamQuat &, const Quat &);

G3TimestreamQuat pow(const G3TimestreamQuat &a, int b);

G3MAP_OF(std::string, G3VectorQuat, G3MapVectorQuat);
G3MAP_OF(std::string, Quat, G3MapQuat);

#endif
