#include <pybindings.h>
#include <container_pybindings.h>
#include <maps/pointing.h>
#include <G3Map.h>
#include <G3Units.h>

#include <boost/math/constants/constants.hpp>

#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <time.h>

using namespace boost::math;

const double PI = constants::pi<double>();

#define ASIN asin
#define ATAN2 atan2

#define QNORM(q) \
	{ \
		double n = dot3(q, q); \
		if (fabs(n - 1.0) > 1e-6) \
			q /= sqrt(n); \
	}

/*
 * Quaternions cannot represent parity flips.  Since celestial coordinates 
 * and az-el coordinates by construction have a different parity, we can't use
 * the general alpha delta angle to x-y-z mapping for one of the coordiate
 * systems.
 *
 * For the Euclidean quaternion representation at the pole,
 * the z coordinate = -sin(elevation) = sin(declination)
 */

static quat
project_on_plane(quat plane_normal, quat point)
{
	// Projects the quaternion onto a plane with unit normal plane_normal
	//   The plane is defined as going through the origin 
	//   with normal = plane_normal

	quat out_q(point);
	//ensure unit vec
	QNORM(plane_normal);
	out_q -= plane_normal * dot3(plane_normal, point);
	QNORM(out_q);
	return out_q;
}

quat
ang_to_quat(double alpha, double delta)
{
	double c_delta = cos(delta / G3Units::rad);
	return quat(0, 
		    c_delta * cos(alpha/G3Units::rad),
		    c_delta * sin(alpha/G3Units::rad),
		    sin(delta / G3Units::rad));
}

void
quat_to_ang(quat q, double &alpha, double &delta)
{
	QNORM(q);
	delta = ASIN(q.R_component_4()) * G3Units::rad;
	alpha = ATAN2(q.R_component_3(), q.R_component_2())*G3Units::rad;
	if (alpha < 0)
		alpha += 360 * G3Units::deg;
}

static boost::python::tuple
py_quat_to_ang(quat q)
{
	double a,d;
	quat_to_ang(q, a, d);

	return boost::python::make_tuple(a, d);
}

double
quat_ang_sep(quat a, quat b)
{
	QNORM(a);
	QNORM(b);

	double d = dot3(a, b);
	if (d > 1)
		return 0;
	if (d < -1)
		return PI * G3Units::rad;
	return acos(d) * G3Units::rad;
}

static quat
coord_quat_to_delta_hat(quat q)
{
	// computes the delta hat vector for a given point on the unit sphere
	// specified by q
	// 
	// (The delta hat is equal to -alpha hat)
	QNORM(q);
	double st = sqrt(1 - (q.R_component_4()*q.R_component_4()));
	quat u= quat(0, 
		     -1 * (q.R_component_2() * q.R_component_4())/st,
		     -1 * (q.R_component_3() * q.R_component_4())/st,
		     st);
	QNORM(u);
	return u;
}

static double
get_rot_ang(quat start_q, quat trans)
{
	// delta is the physicist spherical coordinates delta
	// Computes delta hat for the start q applies trans to it
	// and then computes the angle between that and end_q's delta hat.

	quat t = trans * coord_quat_to_delta_hat(start_q) * ~trans;
	quat end_q = trans * start_q * ~trans;
	quat t_p = coord_quat_to_delta_hat(end_q);
	double sf = (dot3(end_q, cross3(t, t_p)) < 0) ? -1 : 1;
	return sf * quat_ang_sep(t, t_p);
}


static quat
get_transform_quat(double as_0, double ds_0, double ae_0, double de_0,
    double as_1, double ds_1, double ae_1, double de_1)
{
	/*
	 * as = alpha start
	 * ds = delta start
	 * ae = alpha end
	 * de = delta end
	 *
	 * The numeral postscripts are for which set of points.
	 *
	 * Computes a rotation that will take: (as_0,ds_0) to (ae_0, de_0) and
	 * (as_1, ds_1) to (ae_1, de_1)
	 *
	 */

	quat asds_0 = ang_to_quat(as_0, ds_0);
	quat asds_1 = ang_to_quat(as_1, ds_1);
	quat aede_0 = ang_to_quat(ae_0, de_0);
	quat aede_1 = ang_to_quat(ae_1, de_1);

	quat tquat = cross3(asds_0, aede_0);
	double mag = sqrt(dot3(tquat, tquat));
	double ang = quat_ang_sep(asds_0, aede_0);
	tquat *= sin(ang/2.0) / mag;
	tquat += quat(cos(ang/2.0),0,0,0);

	// trans_asds_1 and aede_1 should now be the same up to a rotation
	// around aede_0
	quat trans_asds_1 = tquat * asds_1 * ~tquat;

	// Project them on to a plane and find the angle between the two vectors
	// using (ae_0, de_0) as the normal since we are rotating around that
	// vector.
	quat p_asds1 = project_on_plane(aede_0, trans_asds_1);	
	quat p_aede1 = project_on_plane(aede_0, aede_1);

	double rot_ang = quat_ang_sep(p_asds1, p_aede1);
	double sf = (dot3(aede_0, cross3(p_asds1, p_aede1)) < 0) ? -1 : 1;
	rot_ang *= sf;

	double sin_rot_ang_ov_2 = sin(rot_ang/2.0);
	quat rot_quat = quat(cos(rot_ang/2.0),
			     sin_rot_ang_ov_2 * aede_0.R_component_2(), 
			     sin_rot_ang_ov_2 * aede_0.R_component_3(), 
			     sin_rot_ang_ov_2 * aede_0.R_component_4());
	quat final_trans = rot_quat * tquat;

	return final_trans;
}

quat
get_fk5_j2000_to_gal_quat()
{
	// returns the quaternion that rotates FK5 J2000 to galactic J2000 coordinates
	// return get_transform_quat(0,0, 1.6814025470759737, -1.050488399695429,
	//     0,-0.7853981633974483, 5.750520098164818, -1.2109809382060603);
	return quat(0.4889475076,-0.483210684,0.1962537583,0.699229742);
}

quat
offsets_to_quat(double x_offset, double y_offset)
{
	// Rotates the point (1,0,0) by the rotation matrix for the y_offset
	// and then the rotation matrix for the x_offset
	// quat t = (quat(cos(x_offset/(2.0*G3Units::rad)),0,0,sin(x_offset/(2.0*G3Units::rad))) *
	// 	  quat(cos(y_offset/(2.0*G3Units::rad)),0,sin(y_offset/(2.0*G3Units::rad)),0));
	// return t*quat(0,1.0,0,0)/t;
	// The above is exactly equal to:
	return ang_to_quat(x_offset, -y_offset);
}

quat
get_origin_rotator(double alpha, double delta)
{
	// Rotates the point (1,0,0) to the point specified by alpha and
	// delta via a rotation about the y axis and then the z axis
        return (quat(cos(alpha/2.0), 0, 0, sin(alpha/2.0)) *
                quat(cos(delta/2.0), 0, -sin(delta/2.0), 0));
}

static G3TimestreamQuat
get_origin_rotator_timestream(const G3Timestream &alpha, const G3Timestream &delta,
    MapCoordReference coord_sys)
{
	// Creates the transform that takes (1,0,0) to az, -el 
	// for why it's -el see the comment at the top of this document

	g3_assert(alpha.size() == delta.size());
	G3TimestreamQuat trans_quats(alpha.size(), quat(1,0,0,0));
	trans_quats.start = alpha.start;
	trans_quats.stop = alpha.stop;
	if (coord_sys == Local)
		for (size_t i = 0; i < alpha.size(); i++)
			trans_quats[i] = get_origin_rotator(alpha[i], -delta[i]);
	else
		for (size_t i = 0; i < alpha.size(); i++)
			trans_quats[i] = get_origin_rotator(alpha[i], delta[i]);

	return trans_quats;
}

static G3TimestreamQuat
get_boresight_rotator_timestream(const G3Timestream &az_0, const G3Timestream &el_0,
     const G3Timestream &ra_0, const G3Timestream &dec_0, 
     const G3Timestream &az_1, const G3Timestream &el_1, 
     const G3Timestream &ra_1, const G3Timestream &dec_1)
{
	// Computes the transform that takes (1,0,0) to the point (ra_0, dec_0)
	// and properly handles rotation about the (ra_0, dec_0) point with the
	// inclusion of the second set of points.
	//
	// Stores the output in trans_quats.

	g3_assert(az_0.size() == el_0.size());
	g3_assert(az_0.size() == el_1.size());
	g3_assert(az_0.size() == az_1.size());
	g3_assert(az_0.size() == dec_0.size());
	g3_assert(az_0.size() == dec_1.size());
	g3_assert(az_0.size() == ra_0.size());
	g3_assert(az_0.size() == ra_1.size());
	G3TimestreamQuat trans_quats(ra_0.size(), quat(1,0,0,0));
	trans_quats.start = az_0.start;
	trans_quats.stop = az_0.stop;

	for (size_t i = 0; i < ra_0.size(); i++) {
		trans_quats[i] = get_transform_quat(
		    az_0[i], -el_0[i],
		    ra_0[i], dec_0[i],
		    az_1[i], -el_1[i],
		    ra_1[i], dec_1[i]
		    )*get_origin_rotator(az_0[i], -el_0[i]);
	}

	return trans_quats;
}

G3VectorQuat
get_detector_pointing_quats(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, MapCoordReference coord_sys)
{
	quat q_off = offsets_to_quat(x_offset, y_offset);
	size_t nsamp = trans_quat.size();
	G3VectorQuat det_quats(nsamp, quat(0, 1, 0, 0));

	for (size_t i = 0; i < nsamp; i++)
		det_quats[i] = trans_quat[i] * q_off * ~trans_quat[i];

	if (coord_sys == Local) {
		for (size_t i = 0; i < nsamp; i++) {
			const quat &q = det_quats[i];
			det_quats[i] = quat(q.R_component_1(), q.R_component_2(),
			    q.R_component_3(), -q.R_component_4());
		}
	}

	return det_quats;
}

std::vector<size_t>
get_detector_pointing_pixels(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, G3SkyMapConstPtr skymap)
{
	quat q_off = offsets_to_quat(x_offset, y_offset);
	size_t nsamp = trans_quat.size();
	std::vector<size_t> pixels(nsamp, (size_t) -1);
	quat q;

	if (skymap->coord_ref == Local) {
		for (size_t i = 0; i < nsamp; i++) {
			q = trans_quat[i] * q_off * ~trans_quat[i];
			q = quat(q.R_component_1(), q.R_component_2(),
			    q.R_component_3(), -q.R_component_4());
			pixels[i] = skymap->QuatToPixel(q);
		}
	} else {
		for (size_t i = 0; i < nsamp; i++) {
			q = trans_quat[i] * q_off * ~trans_quat[i];
			pixels[i] = skymap->QuatToPixel(q);
		}
	}

	return pixels;
}

void
get_detector_pointing(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, MapCoordReference coord_sys,
    std::vector<double> &alpha, std::vector<double> &delta)
{
	// For a detector x/y offset and a boresight position specified by
	// trans_quat with a given coordinate system coord_sys,
	// computes the individual detector pointing coordinates.

	quat det_pos = offsets_to_quat(x_offset, y_offset);
	delta.resize(trans_quat.size());
	alpha.resize(trans_quat.size());

	if ((!std::isfinite(x_offset)) || (!std::isfinite(y_offset))){
		log_debug("Found non-finite (inf or nan) offsets");
		for (size_t i=0; i<alpha.size(); i++){
			alpha[i] = nan("");
			delta[i] = nan("");
		}
		return;
	}

	for (size_t i = 0; i < alpha.size(); i++) {
		//uses an inverse that assumes we are on the unit sphere
		quat q = trans_quat[i] * det_pos * ~trans_quat[i];
		quat_to_ang(q, alpha[i], delta[i]);
	}
	if (coord_sys == Local) {
		for (size_t i = 0; i < delta.size(); i++)
			delta[i] *= -1;
	}

}

std::vector<double>
get_detector_rotation(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat)
{
	// Computes the polarization angle rotation that occurs under the
	// transform trans_quat and stores it in rot.
	std::vector<double> rot(trans_quat.size(), 0);
	quat det_pos = offsets_to_quat(x_offset, y_offset);
	for (size_t i = 0; i < rot.size(); i++)
		rot[i] = get_rot_ang(det_pos, trans_quat[i]);

	return rot;
}

PYBINDINGS("maps")
{
	using namespace boost::python;

	// for testing
	def("c_quat_to_ang_", py_quat_to_ang);
	def("c_ang_to_quat_", ang_to_quat);

	def("get_fk5_j2000_to_gal_quat", get_fk5_j2000_to_gal_quat,
	    "Return the rotation quaternion to rotate from equatorial to "
	    "galactic coordinates.");
	def("get_origin_rotator", get_origin_rotator, (arg("alpha"), arg("delta")),
	    "Compute the transformation quaternion that would rotate the "
	    "vector (1, 0, 0) to point in the given direction.");
	def("offsets_to_quat", offsets_to_quat, (arg("x"), arg("y")),
	    "Returns the vector quaternion (0,1,0,0) rotated by the given "
	    "x and y offsets.  Equivalent to ``t * quat(0,1,0,0) / t``, where "
	    "``t = get_origin_rotator(x, -y)``");
	def("get_transform_quat", get_transform_quat,
	    "Computes a rotation that will take (as_0,ds_0) to (ae_0, de_0) and "
	    "(as_1, ds_1) to (ae_1, de_1)");
	def("get_rot_ang", get_rot_ang, (arg("start_q"), arg("trans")),
	    "Computes the boresight rotation of the vector start_q when rotated "
	    "by trans.");
	def("get_origin_rotator_timestream", get_origin_rotator_timestream,
	    (arg("alpha"), arg("delta"), arg("coord_sys")),
	    "Construct a transform quaternion timestream from timestreams of sky "
	    "coordinates. Equivalent to ``R_z(alpha) * R_y(delta)``.");
	def("get_boresight_rotator_timestream", get_boresight_rotator_timestream,
	    "Construct a transform quaternion timestream from timestreams of "
	    "local and equatorial boresight pointing coordinates.  Computes the "
	    "transform from local (az_0, el_0) coordinates to equatorial "
	    "(ra_0, dec_0), accounting for rotation about the boresight by "
	    "including the second set of points.");
}
