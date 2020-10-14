#ifndef _CORE_G3QUAT_H
#define _CORE_G3QUAT_H

#include <G3Vector.h>

#include <boost/math/quaternion.hpp>
#include <cereal/types/vector.hpp>

typedef boost::math::quaternion<double> quat;

namespace cereal
{
// Define cereal serialization for the Quaternions
template<class A>
void serialize(A & ar, quat & q, unsigned version)
{	
	using namespace cereal;
	double a, b, c, d;
	a = q.R_component_1();
	b = q.R_component_2();
	c = q.R_component_3();
	d = q.R_component_4();
	ar & make_nvp("a", a);
	ar & make_nvp("b", b);
	ar & make_nvp("c", c);
	ar & make_nvp("d", d);
	q = quat(a,b,c,d);
}
}

quat cross3(quat a, quat b);
double dot3(quat a, quat b);

G3VECTOR_OF(quat, G3VectorQuat);

class G3TimestreamQuat : public G3VectorQuat
{
public:
	G3TimestreamQuat() : G3VectorQuat() {}
        G3TimestreamQuat(std::vector<quat>::size_type s) : G3VectorQuat(s) {}
        G3TimestreamQuat(std::vector<quat>::size_type s,
            const quat &val) : G3VectorQuat(s, val) {}
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

namespace boost {
namespace math {
	quat operator ~(quat);
};
};

G3VectorQuat operator ~ (const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, double);
G3VectorQuat &operator *= (G3VectorQuat &, double);
G3VectorQuat operator / (const G3VectorQuat &, double);
G3VectorQuat operator / (double, const G3VectorQuat &);
G3VectorQuat operator / (const G3VectorQuat &, const quat &);
G3VectorQuat operator / (const quat &, const G3VectorQuat &);
G3VectorQuat operator / (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator /= (G3VectorQuat &, double);
G3VectorQuat &operator /= (G3VectorQuat &, const quat &);
G3VectorQuat &operator /= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (double, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, quat);
G3VectorQuat operator * (quat, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, quat);

G3VectorQuat pow(const G3VectorQuat &a, double b);
G3VectorQuat pow(const G3VectorQuat &a, int b);

G3TimestreamQuat operator ~ (const G3TimestreamQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, double);
G3TimestreamQuat operator * (double, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, double);
G3TimestreamQuat operator / (double, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, const quat &);
G3TimestreamQuat operator / (const quat &, const G3TimestreamQuat &);
G3TimestreamQuat operator / (const G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, double);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, const quat &);
G3TimestreamQuat &operator /= (G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat &operator *= (G3TimestreamQuat &, const G3VectorQuat &);
G3TimestreamQuat operator * (double, const G3TimestreamQuat &);
G3TimestreamQuat operator * (const G3TimestreamQuat &, quat);
G3TimestreamQuat operator * (quat, const G3TimestreamQuat &);
G3TimestreamQuat &operator *= (G3TimestreamQuat &, quat);

G3TimestreamQuat pow(const G3TimestreamQuat &a, double b);
G3TimestreamQuat pow(const G3TimestreamQuat &a, int b);

#endif
