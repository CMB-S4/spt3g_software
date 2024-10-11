#ifndef _CORE_G3QUAT_H
#define _CORE_G3QUAT_H

#include <G3Vector.h>
#include <G3Map.h>

class G3Quat : public G3FrameObject
{
public:
	G3Quat() : buf_{0, 0, 0, 0}, versor_(false) {}
	G3Quat(double a, double b, double c, double d, bool versor=false) :
	    buf_{a, b, c, d}, versor_(false) {}
	G3Quat(const G3Quat &q) : buf_{q.a(), q.b(), q.c(), q.d()},
	    versor_(q.versor_) {}

	double a() const { return buf_[0]; }
	double b() const { return buf_[1]; }
	double c() const { return buf_[2]; }
	double d() const { return buf_[3]; }

	bool is_versor() const { return versor_; }
	G3Quat versor() const;

	double real() const;
	G3Quat unreal() const;
	G3Quat conj() const;
	double norm() const;
	double abs() const;

	G3Quat operator ~() const;

	G3Quat &operator +=(const G3Quat &);
	G3Quat &operator -=(const G3Quat &);
	G3Quat &operator *=(double);
	G3Quat &operator *=(const G3Quat &);
	G3Quat &operator /=(double);
	G3Quat &operator /=(const G3Quat &);

	G3Quat operator +(const G3Quat &) const;
	G3Quat operator -(const G3Quat &) const;
	G3Quat operator *(double) const;
	G3Quat operator *(const G3Quat &) const;
	G3Quat operator /(double) const;
	G3Quat operator /(const G3Quat &) const;

	bool operator ==(const G3Quat &) const;
	bool operator !=(const G3Quat &) const;

	template <class A> void serialize(A &ar, unsigned v);

	void * buffer();

	std::string Description() const;
	std::string Summary() const { return Description(); }

private:
	double buf_[4];
	bool versor_;

	void versor_inplace();

	SET_LOGGER("G3Quat");
};

namespace cereal {
	template <class A> struct specialize<A, G3Quat, cereal::specialization::member_serialize> {};
}

G3_POINTERS(G3Quat);
G3_SERIALIZABLE(G3Quat, 1);

G3Quat operator *(double, const G3Quat &);
G3Quat operator /(double, const G3Quat &);

inline double real(const G3Quat &q) { return q.real(); };
inline G3Quat unreal(const G3Quat &q) { return q.unreal(); };
inline G3Quat conj(const G3Quat &q) { return q.conj(); };
inline double norm(const G3Quat &q) { return q.norm(); }
inline double abs(const G3Quat &q) { return q.abs(); }

G3Quat pow(const G3Quat &, int);

G3Quat cross3(const G3Quat &a, const G3Quat &b);
double dot3(const G3Quat &a, const G3Quat &b);

G3VECTOR_SPLIT(G3Quat, G3VectorQuat, 2);

class G3TimestreamQuat : public G3VectorQuat
{
public:
	G3TimestreamQuat() : G3VectorQuat() {}
        G3TimestreamQuat(std::vector<G3Quat>::size_type s) : G3VectorQuat(s) {}
        G3TimestreamQuat(std::vector<G3Quat>::size_type s,
            const G3Quat &val) : G3VectorQuat(s, val) {}
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
G3VectorQuat operator / (const G3VectorQuat &, const G3Quat &);
G3VectorQuat operator / (const G3Quat &, const G3VectorQuat &);
G3VectorQuat operator / (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator /= (G3VectorQuat &, double);
G3VectorQuat &operator /= (G3VectorQuat &, const G3Quat &);
G3VectorQuat &operator /= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (double, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, const G3Quat &);
G3VectorQuat operator * (const G3Quat &, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, const G3Quat &);

G3VectorQuat pow(const G3VectorQuat &a, int b);

G3TimestreamQuat operator ~ (const G3TimestreamQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, double);
G3TimestreamQuat operator * (double, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, double);
G3TimestreamQuat operator / (double, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, const G3Quat &);
G3TimestreamQuat operator / (const G3Quat &, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, double);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, const G3Quat &);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat &operator *= (G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat operator * (double, const G3TimestreamQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, const G3Quat &);
G3TimestreamQuat operator * (const G3Quat &, const G3TimestreamQuat &);
G3TimestreamQuat &operator *= (G3TimestreamQuat &, const G3Quat &);

G3TimestreamQuat pow(const G3TimestreamQuat &a, int b);

G3MAP_OF(std::string, G3VectorQuat, G3MapVectorQuat);
G3MAP_SPLIT(std::string, G3Quat, G3MapQuat, 2);

#endif
