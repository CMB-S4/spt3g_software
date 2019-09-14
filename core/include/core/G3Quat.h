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

G3VectorQuat operator * (const G3VectorQuat &, double);
G3VectorQuat &operator *= (G3VectorQuat &, double);
G3VectorQuat operator / (const G3VectorQuat &, double);
G3VectorQuat operator / (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator /= (G3VectorQuat &, double);
G3VectorQuat &operator /= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, const G3VectorQuat &);
G3VectorQuat operator * (double, const G3VectorQuat &);
G3VectorQuat operator * (const G3VectorQuat &, quat);
G3VectorQuat operator * (quat, const G3VectorQuat &);
G3VectorQuat &operator *= (G3VectorQuat &, quat);

G3VectorQuat pow(const G3VectorQuat &a, double b);
G3VectorQuat pow(const G3VectorQuat &a, int b);

#endif
