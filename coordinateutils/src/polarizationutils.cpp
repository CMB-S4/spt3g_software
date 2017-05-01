#include <pybindings.h>
#include <G3Logging.h>
#include <G3Vector.h>

#include <coordinateutils/polarizationutils.h>

namespace bp = boost::python;

void
pol_qu_gal_to_eq(double l, double b, double q_gal, double u_gal,
                    double &q_fk5_out, double &u_fk5_out)
{
	double rot_sin, rot_cos, theta;

	theta = pol_angle_gal_to_eq(l, b, 0.0);
	rot_sin = sin(2.0*theta);
	rot_cos = cos(2.0*theta);

	q_fk5_out = q_gal*rot_cos + u_gal*rot_sin;
	u_fk5_out = u_gal*rot_cos - q_gal*rot_sin;
}

/* Probably not the best exposure. Really just here for unit testing. */
bp::tuple
pol_qu_gal_to_eq_py(double l, double b, double q_gal, double u_gal)
{
	bp::tuple qu_out;
	double q_out, u_out;

	pol_qu_gal_to_eq(l, b, q_gal, u_gal, q_out, u_out);

	qu_out = bp::make_tuple(q_out, u_out);
	return qu_out;
}

double
pol_angle_gal_to_eq(double l, double b, double pol_in)
{
	/*
	 * Computes polarization rotation from Galactic to Equatorial coordinates.
	 */

	/* Galactic coordinates of Equatorial pole */
	static const double eq_l = 122.93191813706228*M_PI/180.0;
	static const double eq_theta = M_PI/2.0 - 27.128251257939215*M_PI/180.0;
	return pol_angle_coord_transform(eq_theta, eq_l, M_PI/2.0-b, l, pol_in);
}

G3VectorDoublePtr
pol_angle_gal_to_eq_vect(std::vector<double> &l, std::vector<double> &b,
                         std::vector<double> &pol_in)
{
	size_t size = pol_in.size();

	G3VectorDoublePtr pol_out = G3VectorDoublePtr(new G3VectorDouble(size));
	pol_angle_gal_to_eq_arr((double *) &l[0], (double *) &b[0],
	                        (double *) &pol_in[0], (double *) &pol_out->front(),
	                        size);

	return pol_out;
}

void
pol_angle_gal_to_eq_arr(double *l, double *b, double *pol_in, double *pol_out,
                        size_t size)
{
	#ifdef OPENMP_FOUND
	#pragma omp parallel for
	#endif
	for (size_t i = 0; i<size; i++)
		pol_out[i] = pol_angle_gal_to_eq(l[i], b[i], pol_in[i]);
}

template <typename T>
bp::list
pol_angle_gal_to_eq_py(T l, T b, T pol_in)
{
	/*
	 * Computes polarization rotation from Galactic to Equatorial coordinates.
	 */
	bp::list pol_out = bp::list();

	size_t i;
	size_t size = bp::len(l);

	std::vector<double> d_b(size);
	std::vector<double> d_l(size);
	std::vector<double> d_pol(size);

	/* The GIL is the worst thing. */
	std::copy(bp::stl_input_iterator<double>(l),
	          bp::stl_input_iterator<double>(),
	          d_l.begin());
	std::copy(bp::stl_input_iterator<double>(b),
	          bp::stl_input_iterator<double>(),
	          d_b.begin());
	std::copy(bp::stl_input_iterator<double>(pol_in),
	          bp::stl_input_iterator<double>(),
	          d_pol.begin());

#ifdef OPENMP_FOUND
	/* Boost Python lists are stupid and don't allow setting values at an index.
	 * So, we use a temp array and then populate the list at the end. */
	double *res_pol = (double*) malloc(size*sizeof(double));
	if (!res_pol){
		log_fatal("malloc failed in pol_angle_gal_to_eq_py");
		exit(1);
	}
	#pragma omp parallel for
	for (i=0; i<size; i++)
		res_pol[i] = pol_angle_gal_to_eq(d_l[i], d_b[i], d_pol[i]);

	for (i=0; i<size; i++)
		pol_out.append(res_pol[i]);
	free(res_pol);
#else
	for (i=0; i<size; i++)
		pol_out.append(pol_angle_gal_to_eq(d_l[i], d_b[i], d_pol[i]));
#endif
	return pol_out;
}

double
pol_angle_eq_to_gal(double alpha, double delta, double pol_in)
{
	/*
	 * Computes polarization rotation from Equatorial to Galactic coordinates.
	 */

	/* Equatorial coordinates of Galactic pole */
	static const double gal_alpha = 192.85948098508993*M_PI/180.0;
	static const double gal_theta = M_PI/2.0 - 27.128251257308968*M_PI/180.0;
	return pol_angle_coord_transform(gal_theta, gal_alpha, M_PI/2.0-delta,
	                                 alpha, pol_in);
}

G3VectorDoublePtr
pol_angle_eq_to_gal_vect(std::vector<double> &alpha, std::vector<double> &delta,
                         std::vector<double> &pol_in)
{
	size_t size = pol_in.size();

	G3VectorDoublePtr pol_out = G3VectorDoublePtr(new G3VectorDouble(size));
	pol_angle_eq_to_gal_arr((double *) &alpha[0], (double *) &delta[0],
	                        (double *) &pol_in[0], (double *) &pol_out->front(),
	                        size);

	return pol_out;
}

void
pol_angle_eq_to_gal_arr(double *alpha, double *delta, double *pol_in,
                        double *pol_out, size_t size)
{
	#ifdef OPENMP_FOUND
	#pragma omp parallel for
	#endif
	for (size_t i = 0; i<size; i++){
		pol_out[i] = pol_angle_eq_to_gal(alpha[i], delta[i], pol_in[i]);
	}
}

template <typename T>
bp::list
pol_angle_eq_to_gal_py(T alpha, T delta, T pol_in)
{
	/*
	 * Computes polarization rotation from Equatorial to Galactic coordinates.
	 */
	bp::list pol_out = bp::list();

	size_t i;
	size_t size = bp::len(alpha);

	std::vector<double> d_delta(size);
	std::vector<double> d_alpha(size);
	std::vector<double> d_pol(size);

	/* The GIL is not the best. */
	std::copy(bp::stl_input_iterator<double>(alpha),
	          bp::stl_input_iterator<double>(),
	          d_alpha.begin());
	std::copy(bp::stl_input_iterator<double>(delta),
	          bp::stl_input_iterator<double>(),
	          d_delta.begin());
	std::copy(bp::stl_input_iterator<double>(pol_in),
	          bp::stl_input_iterator<double>(),
	          d_pol.begin());

#ifdef OPENMP_FOUND
	/* Boost Python lists are stupid and don't allow setting values at an index.
	 * So, we use a temp array and then populate the list at the end. */
	double *res_pol = (double*) malloc(size*sizeof(double));
	if (!res_pol){
		log_fatal("malloc failed in pol_angle_eq_to_gal_py");
		exit(1);
	}
	#pragma omp parallel for
	for (i=0; i<size; i++)
		res_pol[i] = pol_angle_eq_to_gal(d_alpha[i], d_delta[i], d_pol[i]);

	for (i=0; i<size; i++)
		pol_out.append(res_pol[i]);
	free(res_pol);
#else
	for (i=0; i<size; i++)
		pol_out.append(pol_angle_eq_to_gal(d_alpha[i], d_delta[i], d_pol[i]));
#endif
	return pol_out;
}

double
pol_angle_coord_transform(double pole_theta, double pole_phi, double pix_theta,
                          double pix_phi, double pol)
{
	/*
	 * Computes polarization rotation from between a reference coordinate system
	 * and a new coordinate system with pole at position (pole_theta, pole_phi)
	 * in the reference system. Phi is the azimuthal position, and Theta is the
	 * angular distance from the reference pole.
	 */

	double rot_sin, rot_cos, theta, mod_angle;

	/* Full angle can be inferred from the cos and sin. */
	pol_angle_coord_sin_cos(pole_theta, pole_phi, pix_theta, pix_phi, &rot_sin,
	                        &rot_cos);

	if (std::isnan(rot_cos))
		return pol;

	if (rot_sin<0){
		if (rot_cos>M_PI/2.0)
			theta = -M_PI - rot_sin;
		else
			theta = rot_sin;
	} else {
		theta = rot_cos;
	}

	/* The IAU defines polarization angle to increase counterclockwise. */
	mod_angle = fmod(pol-theta, M_PI);
	if (mod_angle<0)
		return mod_angle + M_PI;
	return mod_angle;
}

void
pol_angle_coord_sin_cos(double pole_theta, double pole_phi, double pix_theta,
                        double pix_phi, double *rot_sin, double *rot_cos)
{
	/*
	 * Computes the sin and cos of the rotation angle of polarization under
	 * coordinate transform.
	 */

	double side_ref_pixel, side_new_pixel, side_ref_new, angle_new_pixel,
	       cos_side_new_pixel;

	angle_new_pixel = fmod(pix_phi-pole_phi, 2.0*M_PI);
	side_ref_pixel = pix_theta;
	side_ref_new = pole_theta;
	cos_side_new_pixel = cos(side_ref_new)*cos(side_ref_pixel) +
	                     sin(side_ref_new)*sin(side_ref_pixel)*cos(angle_new_pixel);
	if (angle_new_pixel > M_PI)
		side_new_pixel = 2.0*M_PI - acos(cos_side_new_pixel);
	else
		side_new_pixel = acos(cos_side_new_pixel);
	*rot_sin = asin(sin(side_ref_new)*sin(angle_new_pixel)/sin(side_new_pixel));
	*rot_cos = acos((cos(side_ref_new)-cos(side_ref_pixel)*cos(side_new_pixel))/
	                (sin(side_ref_pixel)*sin(side_new_pixel)));
}

PYBINDINGS("coordinateutils"){
	bp::def("pol_angle_gal_to_eq", pol_angle_gal_to_eq,
	        (bp::arg("l"), bp::arg("b"), bp::arg("pol_in")),
	        "Compute polarization rotation from Galactic to Equatorial Coords");
	bp::def("pol_angle_gal_to_eq", pol_angle_gal_to_eq_py<bp::list&>,
	        (bp::arg("l"), bp::arg("b"), bp::arg("pol_in")),
	        "Compute polarization rotation from Galactic to Equatorial Coords");
	bp::def("pol_angle_gal_to_eq", pol_angle_gal_to_eq_vect,
	        (bp::arg("l"), bp::arg("b"), bp::arg("pol_in")),
	        "Compute polarization rotation from Galactic to Equatorial Coords");
	bp::def("pol_qu_gal_to_eq", pol_qu_gal_to_eq_py,
	        (bp::arg("l"), bp::arg("b"), bp::arg("Q_gal"), bp::arg("U_gal")),
	        "Convert Q/U Stokes parameters from Galactic to Equatorial");
	bp::def("pol_angle_eq_to_gal", pol_angle_eq_to_gal,
	        (bp::arg("alpha"), bp::arg("delta"), bp::arg("pol_in")),
	        "Compute polarization rotation from Equatorial to Galactic Coords");
	bp::def("pol_angle_eq_to_gal", pol_angle_eq_to_gal_py<bp::list&>,
	        (bp::arg("alpha"), bp::arg("delta"), bp::arg("pol_in")),
	        "Compute polarization rotation from Equatorial to Galactic Coords");
	bp::def("pol_angle_eq_to_gal", pol_angle_eq_to_gal_vect,
	        (bp::arg("alpha"), bp::arg("delta"), bp::arg("pol_in")),
	        "Compute polarization rotation from Equatorial to Galactic Coords");
	bp::def("pol_angle_coord_transform", pol_angle_coord_transform,
	        (bp::arg("pole_theta"), bp::arg("pole_phi"), bp::arg("pix_theta"),
	         bp::arg("pix_phi"), bp::arg("pol")));
}
