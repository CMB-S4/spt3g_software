#ifndef _COORDINATEUTILS_QUATERNION_H
#define _COORDINATEUTILS_QUATERNION_H

#include <G3Frame.h>
#include <G3SkyMap.h>
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

// XXX: Define what the following functions do!

void get_detector_pointing( 
	double x_offset, double y_offset,
	const G3VectorQuat & trans_quat,
	MapCoordReference coord_sys,
	std::vector<double> & alpha, std::vector<double> & delta);

void get_detector_rotation( 
	double x_offset, double y_offset,
	const G3VectorQuat & trans_quat,
	std::vector<double> & rot);

quat get_transform_quat(double as_0, double ds_0,
			double ae_0, double de_0,
			double as_1, double ds_1,
			double ae_1, double de_1);
quat offsets_to_quat(double x_offset, double y_offset);
quat get_fk5_j2000_to_gal_quat();
void quat_to_ang(quat q, double & alpha, double & delta);
quat ang_to_quat(double alpha, double delta);
quat coord_quat_to_delta_hat(quat q);
quat get_origin_rotator(double alpha, double delta);
double get_rot_ang(quat start_q, quat end_q, quat trans);

#endif

