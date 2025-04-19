#include <pybindings.h>
#include <container_pybindings.h>
#include <maps/pointing.h>
#include <G3Map.h>
#include <G3Units.h>

#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>

#include <stdlib.h>
#include <time.h>

#define ASIN asin
#define ATAN2 atan2

/*
 * Quaternions cannot represent parity flips.  Since celestial coordinates 
 * and az-el coordinates by construction have a different parity, we can't use
 * the general alpha delta angle to x-y-z mapping for one of the coordiate
 * systems.
 *
 * For the Euclidean quaternion representation at the pole,
 * the z coordinate = -sin(elevation) = sin(declination)
 */

static Quat
unit_vector(const Quat &q)
{
	double n = q.vnorm();
	if (fabs(n - 1.0) > 1e-6)
		return q / sqrt(n);
	return q;
}

static Quat
project_on_plane(const Quat &plane_normal, const Quat &point)
{
	// Projects the quaternion onto a plane with unit normal plane_normal
	//   The plane is defined as going through the origin 
	//   with normal = plane_normal

	Quat out_q(point);
	//ensure unit vec
	auto un = unit_vector(plane_normal);
	out_q -= un * dot3(un, point);
	return unit_vector(out_q);
}

Quat
ang_to_quat(double alpha, double delta)
{
	double c_delta = cos(delta / G3Units::rad);
	return Quat(0,
		    c_delta * cos(alpha/G3Units::rad),
		    c_delta * sin(alpha/G3Units::rad),
		    sin(delta / G3Units::rad));
}

void
quat_to_ang(const Quat &q, double &alpha, double &delta)
{
	auto uq = unit_vector(q);
	delta = ASIN(uq.d()) * G3Units::rad;
	alpha = ATAN2(uq.c(), uq.b()) * G3Units::rad;
	if (alpha < 0)
		alpha += 360 * G3Units::deg;
}

static py::tuple
py_quat_to_ang(const Quat &q)
{
	double a,d;
	quat_to_ang(q, a, d);

	return py::make_tuple(a, d);
}

double
quat_ang_sep(const Quat &a, const Quat &b)
{
	auto ua = unit_vector(a);
	auto ub = unit_vector(b);

	double d = dot3(ua, ub);
	if (d > 1)
		return 0;
	if (d < -1)
		return M_PI * G3Units::rad;
	return acos(d) * G3Units::rad;
}

static Quat
coord_quat_to_delta_hat(const Quat &q)
{
	// computes the delta hat vector for a given point on the unit sphere
	// specified by q
	// 
	// (The delta hat is equal to -alpha hat)
	auto uq = unit_vector(q);
	double st = sqrt(1 - uq.d() * uq.d());
	double ct = -1.0 * uq.d() / st;
	Quat qd(0, uq.b() * ct, uq.c() * ct, st);
	return unit_vector(qd);
}

static double
get_rot_ang(const Quat &start_q, const Quat &trans)
{
	// delta is the physicist spherical coordinates delta
	// Computes delta hat for the start q applies trans to it
	// and then computes the angle between that and end_q's delta hat.

	auto t = trans * coord_quat_to_delta_hat(start_q) * ~trans;
	auto end_q = trans * start_q * ~trans;
	auto t_p = coord_quat_to_delta_hat(end_q);
	double sf = (dot3(end_q, cross3(t, t_p)) < 0) ? -1 : 1;
	return sf * quat_ang_sep(t, t_p);
}


static Quat
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

	auto asds_0 = ang_to_quat(as_0, ds_0);
	auto asds_1 = ang_to_quat(as_1, ds_1);
	auto aede_0 = ang_to_quat(ae_0, de_0);
	auto aede_1 = ang_to_quat(ae_1, de_1);

	auto tquat = cross3(asds_0, aede_0);
	double mag = sqrt(tquat.vnorm());
	double ang = quat_ang_sep(asds_0, aede_0);
	tquat *= sin(ang/2.0) / mag;
	tquat += Quat(cos(ang/2.0), 0, 0, 0);

	// trans_asds_1 and aede_1 should now be the same up to a rotation
	// around aede_0
	auto trans_asds_1 = tquat * asds_1 * ~tquat;

	// Project them on to a plane and find the angle between the two vectors
	// using (ae_0, de_0) as the normal since we are rotating around that
	// vector.
	auto p_asds1 = project_on_plane(aede_0, trans_asds_1);
	auto p_aede1 = project_on_plane(aede_0, aede_1);

	double rot_ang = quat_ang_sep(p_asds1, p_aede1);
	double sf = (dot3(aede_0, cross3(p_asds1, p_aede1)) < 0) ? -1 : 1;
	rot_ang *= sf;

	double sin_rot_ang_ov_2 = sin(rot_ang/2.0);
	Quat rot_quat(cos(rot_ang/2.0),
		     sin_rot_ang_ov_2 * aede_0.b(),
		     sin_rot_ang_ov_2 * aede_0.c(),
		     sin_rot_ang_ov_2 * aede_0.d());
	return rot_quat * tquat;
}

Quat
get_fk5_j2000_to_gal_quat()
{
	// returns the quaternion that rotates FK5 J2000 to galactic J2000 coordinates
	// return get_transform_quat(0,0, 1.6814025470759737, -1.050488399695429,
	//     0,-0.7853981633974483, 5.750520098164818, -1.2109809382060603);
	return Quat(0.4889475076,-0.483210684,0.1962537583,0.699229742);
}

Quat
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

Quat
get_origin_rotator(double alpha, double delta)
{
	// Rotates the point (1,0,0) to the point specified by alpha and
	// delta via a rotation about the y axis and then the z axis
        return (Quat(cos(alpha/2.0), 0, 0, sin(alpha/2.0)) *
                Quat(cos(delta/2.0), 0, -sin(delta/2.0), 0));
}

static G3TimestreamQuat
get_origin_rotator_timestream(const G3Timestream &alpha, const G3Timestream &delta,
    MapCoordReference coord_sys)
{
	// Creates the transform that takes (1,0,0) to az, -el 
	// for why it's -el see the comment at the top of this document

	g3_assert(alpha.size() == delta.size());
	G3TimestreamQuat trans_quats(alpha.size(), Quat(1,0,0,0));
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
	G3TimestreamQuat trans_quats(ra_0.size(), Quat(1,0,0,0));
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
	auto q_off = offsets_to_quat(x_offset, y_offset);
	size_t nsamp = trans_quat.size();
	G3VectorQuat det_quats(nsamp, Quat(0, 1, 0, 0));

	for (size_t i = 0; i < nsamp; i++)
		det_quats[i] = trans_quat[i] * q_off * ~trans_quat[i];

	if (coord_sys == Local) {
		for (size_t i = 0; i < nsamp; i++) {
			const auto &q = det_quats[i];
			det_quats[i] = Quat(q.a(), q.b(), q.c(), -q.d());
		}
	}

	return det_quats;
}

std::vector<uint64_t>
get_detector_pointing_pixels(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, G3SkyMapConstPtr skymap)
{
	auto q_off = offsets_to_quat(x_offset, y_offset);
	size_t nsamp = trans_quat.size();
	std::vector<uint64_t> pixels(nsamp, (uint64_t) -1);
	Quat q;

	if (skymap->coord_ref == Local) {
		for (size_t i = 0; i < nsamp; i++) {
			q = trans_quat[i] * q_off * ~trans_quat[i];
			q = Quat(q.a(), q.b(), q.c(), -q.d());
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

	auto det_pos = offsets_to_quat(x_offset, y_offset);
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
		auto q = trans_quat[i] * det_pos * ~trans_quat[i];
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
	auto det_pos = offsets_to_quat(x_offset, y_offset);
	for (size_t i = 0; i < rot.size(); i++)
		rot[i] = get_rot_ang(det_pos, trans_quat[i]);

	return rot;
}

PYBINDINGS("maps", scope)
{
	// for testing
	scope.def("c_quat_to_ang_", py_quat_to_ang);
	scope.def("c_ang_to_quat_", ang_to_quat);

	scope.def("get_fk5_j2000_to_gal_quat", get_fk5_j2000_to_gal_quat,
	    "Return the rotation quaternion to rotate from equatorial to "
	    "galactic coordinates.");
	scope.def("get_origin_rotator", get_origin_rotator, py::arg("alpha"), py::arg("delta"),
	    "Compute the transformation quaternion that would rotate the "
	    "vector (1, 0, 0) to point in the given direction.");
	scope.def("offsets_to_quat", offsets_to_quat, py::arg("x"), py::arg("y"),
	    "Returns the vector quaternion (0,1,0,0) rotated by the given "
	    "x and y offsets.  Equivalent to ``t * quat(0,1,0,0) / t``, where "
	    "``t = get_origin_rotator(x, -y)``");
	scope.def("get_transform_quat", get_transform_quat,
	    "Computes a rotation that will take (as_0,ds_0) to (ae_0, de_0) and "
	    "(as_1, ds_1) to (ae_1, de_1)");
	scope.def("get_rot_ang", get_rot_ang, py::arg("start_q"), py::arg("trans"),
	    "Computes the boresight rotation of the vector start_q when rotated "
	    "by trans.");
	scope.def("get_origin_rotator_timestream", get_origin_rotator_timestream,
	    py::arg("alpha"), py::arg("delta"), py::arg("coord_sys"),
	    "Construct a transform quaternion timestream from timestreams of sky "
	    "coordinates. Equivalent to ``R_z(alpha) * R_y(delta)``.");
	scope.def("get_boresight_rotator_timestream", get_boresight_rotator_timestream,
	    "Construct a transform quaternion timestream from timestreams of "
	    "local and equatorial boresight pointing coordinates.  Computes the "
	    "transform from local (az_0, el_0) coordinates to equatorial "
	    "(ra_0, dec_0), accounting for rotation about the boresight by "
	    "including the second set of points.");
}
