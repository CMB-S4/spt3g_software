#include <pybindings.h>
#include <G3Logging.h>
#include <G3Vector.h>

#include <coordinateutils/wcs.h>
#include <coordinateutils/coordinateutils.h>

namespace bp = boost::python;

void
coord_gal_to_eq(double *l_in, double *b_in, double *ra_out, double *dec_out,
                size_t size)
{
	double theta_tmp, phi_tmp;

	#ifdef OPENMP_FOUND
	#pragma omp parallel for private(theta_tmp, phi_tmp)
	#endif
	for(size_t i = 0; i<size; i++){
		theta_tmp = l_in[i]*180.0/M_PI;
		phi_tmp = b_in[i]*180.0/M_PI;
		// wcs treats (theta, phi) as equivalent to (l, b), (ra, dec)
		wcscon(WCS_GALACTIC, WCS_J2000, 0.0, 0.0, &theta_tmp, &phi_tmp, 2000.0);
		ra_out[i] = theta_tmp*M_PI/180.0;
		dec_out[i] = phi_tmp*M_PI/180.0;
	}
}

template <typename T>
bp::tuple
coord_gal_to_eq_py(T l_in, T b_in)
{
	bp::list ra_out;
	bp::list dec_out;
	bp::tuple coords_out;

	size_t size = bp::len(l_in);

	double *l = (double*) malloc(size*sizeof(double));
	double *b = (double*) malloc(size*sizeof(double));
	double *ra = (double*) malloc(size*sizeof(double));
	double *dec = (double*) malloc(size*sizeof(double));

	std::copy(bp::stl_input_iterator<double>(l_in),
	          bp::stl_input_iterator<double>(), l);
	std::copy(bp::stl_input_iterator<double>(b_in),
	          bp::stl_input_iterator<double>(), b);

	coord_gal_to_eq(l, b, ra, dec, size);

	for(size_t i=0; i<size; i++){
		ra_out.append(ra[i]);
		dec_out.append(dec[i]);
	}

	/* I'm free! */
	free(l);
	free(b);
	/* I'm freeee! */
	free(ra);
	free(dec);
	/* And I'm waiting... for you... to follow ME! */

	coords_out = bp::make_tuple(ra_out, dec_out);
	return coords_out;
}

bp::tuple
coord_gal_to_eq_vect(std::vector<double> &l_in, std::vector<double> &b_in)
{
	g3_assert(l_in.size()==b_in.size());
	size_t size = l_in.size();
	bp::tuple coords_out;

	G3VectorDoublePtr ra_out = G3VectorDoublePtr(new G3VectorDouble(size));
	G3VectorDoublePtr dec_out = G3VectorDoublePtr(new G3VectorDouble(size));

	coord_gal_to_eq(&l_in[0], &b_in[0], &ra_out->front(), &dec_out->front(), size);

	coords_out = bp::make_tuple(ra_out, dec_out);
	return coords_out;
}

void
coord_eq_to_gal(double *ra_in, double *dec_in, double *l_out, double *b_out,
                size_t size)
{
	double theta_tmp, phi_tmp;

	#ifdef OPENMP_FOUND
	#pragma omp parallel for private(theta_tmp, phi_tmp)
	#endif
	for(size_t i = 0; i<size; i++){
		theta_tmp = ra_in[i]*180.0/M_PI;
		phi_tmp = dec_in[i]*180.0/M_PI;
		// wcs treats (theta, phi) as equivalent to (l, b), (ra, dec)
		wcscon(WCS_J2000, WCS_GALACTIC, 0.0, 0.0, &theta_tmp, &phi_tmp, 2000.0);
		l_out[i] = theta_tmp*M_PI/180.0;
		b_out[i] = phi_tmp*M_PI/180.0;
	}
}

template <typename T>
bp::tuple
coord_eq_to_gal_py(T ra_in, T dec_in)
{
	bp::list l_out;
	bp::list b_out;
	bp::tuple coords_out;

	size_t size = bp::len(ra_in);

	double *ra = (double*) malloc(size*sizeof(double));
	double *dec = (double*) malloc(size*sizeof(double));
	double *l = (double*) malloc(size*sizeof(double));
	double *b = (double*) malloc(size*sizeof(double));

	/* The GIL sucks. */
	for(size_t i=0; i<size; i++){
		ra[i] = bp::extract<double>(ra_in[i]);
		dec[i] = bp::extract<double>(dec_in[i]);
	}

	coord_eq_to_gal(ra, dec, l, b, size);

	for(size_t i=0; i<size; i++){
		l_out.append(l[i]);
		b_out.append(b[i]);
	}

	/* I'm free! */
	free(l);
	free(b);
	/* I'm freeee! */
	free(ra);
	free(dec);
	/* And I'm waiting... for you... to follow ME! */

	coords_out = bp::make_tuple(l_out, b_out);
	return coords_out;
}

bp::tuple
coord_eq_to_gal_vect(std::vector<double> &ra_in, std::vector<double> &dec_in)
{
	g3_assert(ra_in.size()==dec_in.size());
	size_t size = ra_in.size();
	bp::tuple coords_out;

	G3VectorDoublePtr l_out = G3VectorDoublePtr(new G3VectorDouble(size));
	G3VectorDoublePtr b_out = G3VectorDoublePtr(new G3VectorDouble(size));

	coord_gal_to_eq(&ra_in[0], &dec_in[0], &l_out->front(), &b_out->front(), size);

	coords_out = bp::make_tuple(l_out, b_out);
	return coords_out;
}

PYBINDINGS("coordinateutils"){
	bp::def("coord_gal_to_eq", coord_gal_to_eq_py<bp::list&>,
	        (bp::arg("l"), bp::arg("b")),
	        "Returns (ra, dec) tuple corresponding to input (l, b)");
	bp::def("coord_gal_to_eq", coord_gal_to_eq_vect,
	        (bp::arg("l"), bp::arg("b")),
	        "Returns (ra, dec) tuple corresponding to input (l, b)");
	bp::def("coord_eq_to_gal", coord_eq_to_gal_py<bp::list&>,
	        (bp::arg("ra"), bp::arg("dec")),
	        "Returns (l, b) tuple corresponding to input (ra, dec)");
	bp::def("coord_eq_to_gal", coord_eq_to_gal_vect,
	        (bp::arg("ra"), bp::arg("dec")),
	        "Returns (l, b) tuple corresponding to input (ra, dec)");
}
