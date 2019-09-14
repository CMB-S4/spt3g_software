#include <pybindings.h>
#include <container_pybindings.h>
#include <coordinateutils/pointing.h>
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

//#define CHECK_QUAT_INVERSE


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
	plane_normal /= sqrt(dot3(plane_normal,plane_normal));
	out_q -= plane_normal * dot3(plane_normal, point);
	return out_q;
}

static bool
sloppy_eq(quat a, quat b, double slop = 1e-6)
{
	// Fuzzy quaternion equality comparison
	return ((fabs(a.R_component_1() - (b.R_component_1())) < slop ) &&
		(fabs(a.R_component_2() - (b.R_component_2())) < slop ) &&
		(fabs(a.R_component_3() - (b.R_component_3())) < slop ) &&
		(fabs(a.R_component_4() - (b.R_component_4())) < slop ));
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
	double d = dot3(q,q);
	if (fabs(d - 1.0) > 1e-6){
		q /= sqrt(d);
	}
	delta = ASIN(q.R_component_4()) * G3Units::rad;
	alpha = ATAN2(q.R_component_3(), q.R_component_2())*G3Units::rad;
}

static boost::python::tuple
py_quat_to_ang(quat q)
{
	double a,d;
	quat_to_ang(q, a, d);

	return boost::python::make_tuple(a, d);
}

quat
coord_quat_to_delta_hat(quat q)
{
	// computes the delta hat vector for a given point on the unit sphere
	// specified by q
	// 
	// (The delta hat is equal to -alpha hat)

	q /= sqrt(dot3(q,q));
	double st = sqrt(1 - (q.R_component_4()*q.R_component_4()));
	quat u= quat(0, 
		     -1 * (q.R_component_2() * q.R_component_4())/st,
		     -1 * (q.R_component_3() * q.R_component_4())/st,
		     st);
	u /= sqrt(dot3(u,u));
	return u;
}

double
get_rot_ang(quat start_q, quat end_q, quat trans)
{
	// delta is the physicist spherical coordinates delta
	// Computes delta hat for the start q applies trans to it
	// and then computes the angle between that and end_q's delta hat.

	quat t = trans * coord_quat_to_delta_hat(start_q) / trans;
	quat t_p = coord_quat_to_delta_hat(end_q);
	
	t /= sqrt(dot3(t,t));
	t_p /= sqrt(dot3(t_p,t_p));
	double d  = dot3(t,t_p);
	double sf = (dot3(end_q, cross3(t, t_p)) < 0) ? -1 : 1;
	if (d > 1) {
		g3_assert(d < 1.01);
		return 0;
	} else if (d < -1) {
		g3_assert(d > -1.01);
		return PI * G3Units::rad;
	} else {
		return sf * acos(d) * G3Units::rad;
	}
}


quat
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
	double ang = acos(dot3(asds_0, aede_0)); 
	tquat *= sin(ang/2.0) / mag;
	tquat += quat(cos(ang/2.0),0,0,0);

	// trans_asds_1 and aede_1 should now be the same up to a rotation
	// around aede_0
	quat trans_asds_1 = tquat * asds_1 / tquat;

	// Project them on to a plane and find the angle between the two vectors
	// using (ae_0, de_0) as the normal since we are rotating around that
	// vector.
	quat p_asds1 = project_on_plane(aede_0, trans_asds_1);	
	quat p_aede1 = project_on_plane(aede_0, aede_1);
	p_asds1 /= sqrt(dot3(p_asds1,p_asds1));
	p_aede1 /= sqrt(dot3(p_aede1,p_aede1));

	double rot_ang = acos(dot3(p_asds1, p_aede1));
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

static std::vector<double>
test_trans(double az_0, double el_0, double ra_0, double dec_0,
    double az_1, double el_1, double ra_1, double dec_1,
    double az_t, double el_t)
{
	// computes the transform from the first 4 variables
	// returns that transform to the last 2 variables
	
	double ra_t, dec_t;
	quat q = get_transform_quat(az_0, -el_0,
				    ra_0, dec_0,
				    az_1, -el_1,
				    ra_1, dec_1);
	quat azel = ang_to_quat(az_t, -el_t);
	quat rad = q * azel / q;
	quat_to_ang(rad, ra_t, dec_t);
	std::vector<double> r(2,0);
	r[0] = ra_t;
	r[1] = dec_t;
	return r;
}

static std::vector<double>
test_gal_trans(double az_0, double el_0, double ra_0, double dec_0,
    double az_1, double el_1, double ra_1, double dec_1,
    double az_t, double el_t)
{
	double l_t, b_t;
	quat q = get_fk5_j2000_to_gal_quat() *get_transform_quat(az_0, -el_0,
								 ra_0, dec_0,
								 az_1, -el_1,
								 ra_1, dec_1);
	quat azel = ang_to_quat(az_t, -el_t);
	quat rad = q * azel / q;
	quat_to_ang(rad, l_t, b_t);
	std::vector<double> r(2,0);
	r[0] = l_t;
	r[1] = b_t;
	return r;
}

static double
test_gal_trans_rot(double ra, double dec)
{
	quat start = ang_to_quat(ra,dec);
	quat trans = get_fk5_j2000_to_gal_quat();

	return get_rot_ang(start, trans*start/trans, trans);
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

static void
print_fk5_j2000_to_gal_quat()
{
	// uhh, so, this code was a super lazy way to get the quaternion
	// that takes fk5 j2000 to galactic j2000
	std::cout << std::setprecision(10) << std::endl;
	std::cout << get_transform_quat(
		0,0, 1.6814025470759737, -1.050488399695429,
		0,-0.7853981633974483, 5.750520098164818, -1.2109809382060603)
	          << std::endl;
}

quat
get_fk5_j2000_to_gal_quat()
{
	// returns the quaternion that rotates fk5j2000 to galactic J2000
	// coordinates
	return quat(0.4889475076,-0.483210684,0.1962537583,0.699229742);
}

static void
create_det_az_el_trans(const G3Timestream &az, const G3Timestream &el,
    G3VectorQuat &trans_quats)
{
	// Creates the transform that takes (1,0,0) to az, -el 
	// for why it's -el see the comment at the top of this document

	g3_assert(az.size() == el.size());
	trans_quats = G3VectorQuat(az.size(), quat(1,0,0,0));
	for (size_t i = 0; i < az.size(); i++)
		trans_quats[i] = get_origin_rotator(az[i], -el[i]);
}

static void
create_lazy_det_ra_dec_trans(const G3Timestream &ra, const G3Timestream &dec, 
    G3VectorQuat &trans_quats)
{
	// Creates the transform that takes (1,0,0) to ra,dec
	g3_assert(ra.size() == dec.size());
	trans_quats = G3VectorQuat(ra.size(), quat(1,0,0,0));
	for (size_t i = 0; i < ra.size(); i++)
		trans_quats[i] = get_origin_rotator(ra[i], dec[i]);
}

static void
create_det_ra_dec_trans(const G3Timestream &az_0, const G3Timestream &el_0, 
     const G3Timestream &ra_0, const G3Timestream &dec_0, 
     const G3Timestream &az_1, const G3Timestream &el_1, 
     const G3Timestream &ra_1, const G3Timestream &dec_1, 
     G3VectorQuat & trans_quats)
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
	trans_quats = G3VectorQuat(ra_0.size(), quat(1,0,0,0));	

	for (size_t i = 0; i < ra_0.size(); i++) {
		trans_quats[i] = get_transform_quat(
		    az_0[i], -el_0[i],
		    ra_0[i], dec_0[i],
		    az_1[i], -el_1[i],
		    ra_1[i], dec_1[i]
		    )*get_origin_rotator(az_0[i], -el_0[i]);
	}
}

static void
convert_ra_dec_trans_to_gal(const G3VectorQuat &radec_trans,
    G3VectorQuat &gal_trans)
{
	// Converts a rotation from (1,0,0) to fk5 j2000 into a rotation that
	// takes (1,0,0) to our galactic (l,b)

	gal_trans = G3VectorQuat(radec_trans.size(), quat(1,0,0,0));
	quat gt = get_fk5_j2000_to_gal_quat();
	for (size_t i = 0; i < radec_trans.size(); i++)
		gal_trans[i] = gt*radec_trans[i];
}

void
get_detector_pointing(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, MapCoordReference coord_sys,
    std::vector<double> &alpha, std::vector<double> &delta)
{
	// For a detector x/y offset and a boresight position specified by
	// trans_quat with a given coordinate system coord_sys,
	// computes the individual detector pointing coordinates.
	//
	// Assumes alpha and delta have already been allocated.

	quat det_pos = offsets_to_quat(x_offset, y_offset);

	if ((!std::isfinite(x_offset)) || (!std::isfinite(y_offset))){
		log_debug("Found non-finite (inf or nan) offsets");
		for (size_t i=0; i<alpha.size(); i++){
			alpha[i] = nan("");
			delta[i] = nan("");
		}
		return;
	}

	for (size_t i = 0; i < alpha.size(); i++) {
		//using boost inverse
		//quat q=trans_quat[i]*det_pos/trans_quat[i];
		
		//uses an inverse that assumes we are on the unit sphere
		const quat & t = trans_quat[i];
		quat q=trans_quat[i]*det_pos * quat( t.R_component_1(),
		    -t.R_component_2(), -t.R_component_3(), -t.R_component_4());

		quat_to_ang(q, alpha[i], delta[i]);

		#ifdef CHECK_QUAT_INVERSE
		double a,d;
		quat u = trans_quat[i]*det_pos/trans_quat[i];		
		quat_to_ang(u, a, d);
		if( fabs(a - alpha[i]) > 1e-5 || fabs(d - delta[i]) > 1e-5){
			log_fatal("Failed trans %lf %lf %lf %lf\n", a, alpha[i], d, delta[i]);
		}
		#endif
	}
	if (coord_sys == Local) {
		for (size_t i = 0; i < delta.size(); i++)
			delta[i] *= -1;
	}

}

void
get_detector_rotation(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, std::vector<double> &rot)
{
	// Computes the polarization angle rotation that occurs under the 
	// transform trans_quat and stores it in rot.	   
	
	rot = std::vector<double>(trans_quat.size(), 0);
	quat det_pos = offsets_to_quat(x_offset, y_offset);
	for (size_t i = 0; i < rot.size(); i++) {
		quat q = trans_quat[i]*det_pos/trans_quat[i];
		rot[i] = get_rot_ang(det_pos, q, trans_quat[i]);
	}
}

static std::vector<double>
convert_celestial_offsets_to_local_offsets(G3VectorQuat trans_vec,
    double x_offset_celest, double y_offset_celest)
{
	double alpha, delta;
 	double x_offset_local, y_offset_local;

	// written at a time when we only have python bindings for the vector quat
	g3_assert( trans_vec.size() == 1); 
	quat trans = trans_vec[0];

	// First we compute the celestial position of boresight.
	quat bs_cel = trans * quat(0,1,0,0) / trans;
	quat_to_ang(bs_cel, alpha, delta); //alpha, delta are celestial boresight
	
	// Next, with our celestial position offsets we compute the position
	// of our detector in celestial coordinates with the celestial offsets.
	// because the offset is in celestial coordinates we do that with the transform
	// that takes (1,0,0) in celestial coordinates to the boresight  while ignoring
	//boresight rotation

	quat celestial_trans = get_origin_rotator(alpha, delta);

	quat det_cel = celestial_trans * offsets_to_quat(x_offset_celest, y_offset_celest) / 
	    celestial_trans;

	// We then transform those positions to local coordinates using the full transform
	quat inv_trans = boost::math::conj(trans);
	quat det_local = inv_trans * det_cel / inv_trans;

	#ifdef CHECK_QUAT_INVERSE
	//just to check I'm not insane can be dropped in the future
	g3_assert(sloppy_eq(inv_trans * bs_cel / inv_trans, quat(0,1,0,0)));
	#endif
	
	// we read off the position of the detector.
	quat_to_ang(det_local, x_offset_local, y_offset_local);

	// C++11 is magic.
	return {x_offset_local, -y_offset_local};
}

static double
angle_d2(double ra_0, double dec_0, double ra_1, double dec_1)
{
	double delta_ra_abs = fabs(ra_0 - ra_1);
	double delta_dec = dec_0 - dec_1;
	if (delta_ra_abs > PI){
		delta_ra_abs = 2 * PI - delta_ra_abs;
	}
	delta_ra_abs /= cos(dec_0);
	return delta_ra_abs*delta_ra_abs + delta_dec*delta_dec;

}

static G3VectorQuat
get_closest_transform(double ra, double dec,
    const G3VectorDouble &ras, const G3VectorDouble &decs,
    const G3VectorQuat &trans)
{
	// Returns quaternion transformation nearest the requested ra/dec

	double dist = angle_d2(ras[0]/G3Units::rad, decs[0]/G3Units::rad,
			       ra/G3Units::rad, dec/G3Units::rad);
	size_t ind = 0;
	for (size_t i = 0; i < ras.size(); i++) {
		double td = angle_d2(ras[i]/G3Units::rad, decs[i]/G3Units::rad,
				     ra/G3Units::rad, dec/G3Units::rad);
		if (td < dist) {
			dist = td;
			ind = i;
		}
	}
	G3VectorQuat v;
	v.push_back(trans[ind]);
	return v;
}


PYBINDINGS("coordinateutils")
{
	using namespace boost::python;

	// XXX: All of these need doc strings
	def("test_trans", test_trans);
	def("test_gal_trans", test_gal_trans);
	def("test_gal_trans_rot", test_gal_trans_rot);
	def("print_fk5_j2000_to_gal_quat", print_fk5_j2000_to_gal_quat);

	def("create_det_az_el_trans", create_det_az_el_trans);
	def("create_lazy_det_ra_dec_trans", create_lazy_det_ra_dec_trans);
	def("create_det_ra_dec_trans", create_det_ra_dec_trans);
	def("convert_ra_dec_trans_to_gal", convert_ra_dec_trans_to_gal);
	def("c_quat_to_ang", py_quat_to_ang);
	def("c_ang_to_quat", ang_to_quat);
	def("get_detector_pointing", get_detector_pointing);
	def("get_detector_rotation", get_detector_rotation);
	def("get_origin_rotator", get_origin_rotator);

	def("convert_celestial_offsets_to_local_offsets", 
	    convert_celestial_offsets_to_local_offsets);

	def("get_closest_transform", get_closest_transform);
}

