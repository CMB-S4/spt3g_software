#include <pybindings.h>

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <G3Logging.h>
#include <G3Units.h>
#include <G3SkyMap.h>

#include <coordinateutils/flatskyprojection.h>

#ifdef OPENMP_FOUND
#include <omp.h>
#endif

#define COS cosf
#define SIN sinf
#define ASIN asinf
#define ATAN2 atan2f

using namespace G3Units;

int pixel_2d_to_pixel_1d( int x, int y, int nx, int ny){
	return (x < 0 || y < 0 || x >= nx || y >= ny)  ? nx * ny : x + y * nx;
}

std::vector<int> pixel_1d_to_pixel_2d( int map_index,  int nx, int ny)
{
	std::vector<int> out(2);

	if (map_index >= nx * ny || map_index < 0){
		out[1] = -1;
		out[0] = -1;
	} else {
		out[1] = map_index % nx;
		out[0] = map_index / nx;
	}
	return out;
}

void radec_to_map2d_proj0(double ra,double dec, double ra0, double dec0,
    double res, double x_min, double y_min, int n_x, int n_y, double & xpix, 
    double & ypix) 
{
	double x,y;
	ra = ra - 360*deg * (ra - ra0 > 180*deg) + 
	    360*deg * ( ra - ra0 < -180*deg) ;
	
	x=(ra-ra0) * COS(dec/rad);
	y=(-1.*dec+dec0);
	
	x *= -1;
	y *= -1;
	
	xpix = ((x-x_min)/res);
	ypix = ((y-y_min)/res);
};


void radec_to_map2d_proj1(double ra,double dec, double ra0, double dec0, 
    double res, double x_min, double y_min,  int n_x, int n_y, double & xpix, 
    double & ypix) 
{
	double x,y;
	ra = ra - 360*deg * (ra - ra0 > 180*deg) + 
	    360*deg * ( ra - ra0 < -180*deg) ;
	x=(ra-ra0);
	y=(-1.*dec+dec0);
	
	x *= -1;
	y *= -1;
	
	xpix = ((x-x_min)/res);
	ypix = ((y-y_min)/res);
};


void radec_to_map2d_proj2(double ra,double dec, double ra0, double dec0,
    double res,  double x_min, double y_min, int n_x, int n_y, 
    double negcosdelta0, double sindelta0, double & xpix, double & ypix) 
{
	double x,y;
	double alpha,delta;
	alpha = ra/rad;
	delta = dec/rad;
	double delta0 = (90*deg + dec0)/rad;
	double alpha0 = ra0/rad;

	y = (COS(delta) * negcosdelta0 * COS(alpha - alpha0) - 
	    sindelta0 *SIN(delta)) * rad;
	x = (COS(delta)*SIN(alpha - alpha0)) * rad;

	x *= -1;
	y *= -1;

	xpix = ((x-x_min)/res);
	ypix = ((y-y_min)/res);
}


void radec_to_map2d_proj4(double ra,double dec,  double ra0, double dec0,
    double res, double x_min, double y_min, int n_x, int n_y,
    double negcosdelta0, double sindelta0, double & xpix, double & ypix) 
{
	double x,y;
	double alpha,delta, alpha0, delta0;
	alpha = ra/rad;
	delta = dec/rad;

	alpha0 = ra0/rad;
	delta0 = (90*deg + dec0)/rad;

	double k=rad*(2./(1.+negcosdelta0*SIN(delta) + 
	     sindelta0*COS(delta)*COS(alpha-alpha0)));
	y= k * (COS(delta)*negcosdelta0 *COS(alpha-alpha0) - 
	     sindelta0 *SIN(delta));
	x= k * COS(delta)*SIN(alpha-alpha0);

	x *= -1;
	y *= -1;

	xpix = ((x-x_min)/res);
	ypix = ((y-y_min)/res);
}


void radec_to_map2d_proj5(double ra,double dec, double ra0, double dec0,
    double res, double x_min, double y_min, int n_x, int n_y, 
    double negcosdelta0, double sindelta0, double & xpix, double & ypix) 
{
	double x,y;
	double alpha,delta, alpha0, delta0;
	alpha = ra/rad;
	delta = (dec)/rad;
	alpha0 = ra0/rad;
	delta0 = (90*deg + dec0)/rad;

	double cos_delta = COS(delta);
	double sin_delta = SIN(delta);
	double cos_alphaalpha0 = COS(alpha-alpha0);

	double k =  rad * sqrt(2./(1. + negcosdelta0*sin_delta + 
	    sindelta0*cos_delta*cos_alphaalpha0 ));
	y = -1 * k * (cos_delta * negcosdelta0 * cos_alphaalpha0 - 
	    sindelta0 * sin_delta);
	x = -k * (cos_delta * SIN(alpha-alpha0));
	xpix =  ((x-x_min)/res);
	ypix =  ((y-y_min)/res);
}


int pixel_2d_double_to_1d_index(double xind, double yind, int nx, int ny)
{
	int x =(int) xind;
	int y =(int) yind;
	return (x < 0 || y < 0 || x >= nx || y >= ny)  ? nx * ny : x + y * nx;
}


void angle_to_pixel_1d_vec( const std::vector<double> & ras, 
			    const std::vector<double> & decs,
                            double ra0, double dec0,
                            int n_x, int n_y,
                            double pixel_res, double x_pixel_res,
                            MapProjection proj,
                            std::vector<int> & output_index)
{
	g3_assert(decs.size() == ras.size());
	if (output_index.size() != ras.size()){
		output_index = std::vector<int>(decs.size(), -1);
	}
	double negcosdelta0 =  -1.0 *  COS( (90*deg + dec0)/rad );
	double sindelta0 = SIN((90*deg + dec0)/rad );

	int n_pnts = ras.size();

	if (x_pixel_res <= 0) x_pixel_res = pixel_res;
	if (proj != MapProjection::Proj9) g3_assert(x_pixel_res == pixel_res);

	double res   =   pixel_res;
	double x_res = x_pixel_res;

	double x_min = -0.5*(n_x * res);
	double y_min = -0.5*(n_y * res);

	switch (proj) {
	case Proj0: {
		for (int i=0; i < n_pnts; i++){
			double xpix, ypix;
			radec_to_map2d_proj0(ras[i], decs[i],
					     ra0, dec0, res,
					     x_min, y_min, n_x, n_y,
					     xpix, ypix);
			output_index[i] = 
			    pixel_2d_double_to_1d_index(xpix, ypix, n_x, n_y);
		}
		break;
	}
		
	case Proj1: {
		for (int i=0; i < n_pnts; i++){
			double xpix, ypix;
			radec_to_map2d_proj1(ras[i], decs[i],
			                     ra0, dec0, res,
			                     x_min, y_min, n_x, n_y,
			                     xpix, ypix);
			output_index[i] = 
			    pixel_2d_double_to_1d_index(xpix, ypix, n_x, n_y);
		}
		break;
	}
	case Proj2: {
		for (int i=0; i < n_pnts; i++){
			double xpix, ypix;
			radec_to_map2d_proj2(ras[i], decs[i],
			                     ra0, dec0, res,
			                     x_min, y_min, n_x, n_y,
			                     negcosdelta0, sindelta0,
			                     xpix, ypix );
			output_index[i] = 
			    pixel_2d_double_to_1d_index(xpix, ypix, n_x, n_y);
		}
		break;
	}

	case Proj4: {
		for (int i=0; i < n_pnts; i++){
			double xpix, ypix;
			radec_to_map2d_proj4(ras[i], decs[i],
			                     ra0, dec0, res,
			                     x_min, y_min, n_x, n_y,
			                     negcosdelta0, sindelta0,
			                     xpix, ypix );
			output_index[i] = 
			    pixel_2d_double_to_1d_index(xpix, ypix, n_x, n_y);
		}
		break;
	}

	case Proj5: {
		for (int i=0; i < n_pnts; i++){
			double xpix, ypix;
			radec_to_map2d_proj5(ras[i], decs[i],
			                     ra0, dec0, res,
			                     x_min, y_min, n_x, n_y,
			                     negcosdelta0, sindelta0,
			                     xpix, ypix );
			output_index[i] = 
			    pixel_2d_double_to_1d_index(xpix, ypix, n_x, n_y);
		}
		break;
	}
	default: {
		log_fatal("Err: Proj %d unimplemented...ran out of paste\n",proj);
		break;
	}
	}
}

void pixel_1d_to_angle(const std::vector<int> & pix_ind, double ra0, 
    double dec0, int n_x, int n_y, double pixel_res, double x_res,
    MapProjection proj, bool wrap_ra, std::vector<double> & ra_out, 
    std::vector<double> & dec_out )
{
	size_t npix = pix_ind.size();
	ra_out = std::vector<double>(npix);
	dec_out = std::vector<double>(npix);
	std::vector<double> y_coord(npix);
	std::vector<double> x_coord(npix);

	if (x_res <= 0) x_res = pixel_res;
	if (proj != MapProjection::Proj9) g3_assert(x_res == pixel_res);


	for (size_t i=0; i < npix; i++){
		x_coord[i] = (0.5 * n_x - ((double) (pix_ind[i] % n_x)) - 0.5)
		    * pixel_res/rad;
		y_coord[i] = (0.5 * n_y - ((double) (pix_ind[i] / n_x)) - 0.5) 
		    * pixel_res/rad;
	}

	switch (proj){
	case Proj0:{
		for (size_t i=0; i < npix; i++){
			dec_out[i] = (dec0 - y_coord[i]*rad);
			ra_out[i] = x_coord[i] * rad / COS(dec_out[i]/rad)+ra0;
		}
		break;
	}
	case Proj1:{
		for (size_t i=0; i < npix; i++){
			dec_out[i] = dec0 - y_coord[i]*rad;
			ra_out[i] = x_coord[i]*rad + ra0;
		}
		break;
	}
	case Proj2:{
		for (size_t i=0; i < npix; i++){
			double rho = sqrt( x_coord[i]*x_coord[i] + 
			    y_coord[i]*y_coord[i]);
			double c = ASIN(rho) * rad;
			double alpha_temp = ASIN( COS(c/rad)*SIN(dec0/rad) - 
			    y_coord[i]*COS(dec0/rad)  ) * rad;
			double lambda_temp = 
			    ATAN2(x_coord[i]*SIN(c/rad), 
			    rho*COS(dec0/rad)*COS(c/rad) + 
			    y_coord[i]*SIN(dec0/rad)*SIN(c/rad) ) * rad;
			if (rho < 1e-8){
				alpha_temp = dec0;
				lambda_temp = 0;
			}
			dec_out[i] = alpha_temp;
			ra_out[i] = ra0 + lambda_temp;
		}
		break;
	}
	case Proj4:{
		for (size_t i=0; i < npix; i++){
			double rho = sqrt( x_coord[i]*x_coord[i] + 
			    y_coord[i]*y_coord[i]);
			double c = 2.0 * atan(rho/2.0) * rad;
			double alpha_temp = ASIN(COS(c/rad)*SIN(dec0/rad) - 
			    y_coord[i]*SIN(c/rad)/rho*COS(dec0/rad)) * rad;
			double lambda_temp = rad * ATAN2(x_coord[i]*SIN(c/rad),
			    rho*COS(dec0/rad)*COS(c/rad) + 
			    y_coord[i]*SIN(dec0/rad)*SIN(c/rad));
			if( rho < 1e-8){
				alpha_temp = dec0;
				lambda_temp = 0.0;
			}
			dec_out[i] = alpha_temp;
			ra_out[i] = ra0 + lambda_temp;
		}
		break;
	}
	case Proj5:{
		for (size_t i=0; i < npix; i++){
			double rho = sqrt( x_coord[i]*x_coord[i] + 
			    y_coord[i]*y_coord[i]);
			double c = 2.0*ASIN(rho/2.0) * rad;
			double alpha_temp = rad*ASIN(COS(c/rad)*SIN(dec0/rad) 
			    - y_coord[i]*SIN(c/rad)/rho*COS(dec0/rad)) ;
			double lambda_temp = rad * ATAN2(x_coord[i]*SIN(c/rad),
			    rho*COS(dec0/rad)*COS(c/rad) + 
			    y_coord[i]*SIN(dec0/rad)*SIN(c/rad));
			if( rho < 1e-8){
				alpha_temp = dec0;
				lambda_temp = 0.0;
			}
			dec_out[i] = alpha_temp;
			ra_out[i] = ra0 + lambda_temp;
		}
		break;
	}
	case Proj9:{
		log_fatal("Untested and doesn't account for x resolution "
			  "properly."
			  "You can fix it Jason G.");
		for (size_t i=0; i < npix; i++){
			double pixsize = pixel_res/deg;
			double asx = pixsize/COS( dec0/rad );
			double sx  = round(asx*n_x)/n_x;
			dec_out[i] = dec0 - y_coord[i]*rad;
			ra_out[i]  = x_coord[i]*rad*sx/pixsize + ra0;
		}
		break;
	}
	default: {
		log_fatal("Err: Proj %d unimplemented\n",proj);
		break;
	}
	}

	if (wrap_ra){
		for (size_t i=0; i < npix; i++){
			ra_out[i] = ra_out[i] < 0 ? fmod(360.0*deg * 
			    (1.0 + ceilf(fabs(ra_out[i])/(360.0*deg))) + 
			     ra_out[i], 360.*deg)
			    : fmod(ra_out[i], 360.0*deg);
		}
	} else{
		for (size_t i=0; i < npix; i++){
			ra_out[i] = ra_out[i] - ra0 < -180*deg ? 
			    ra_out[i] + 360.0*deg : ra_out[i];
			ra_out[i] = ra_out[i] - ra0 >= 180*deg ? 
			    ra_out[i] - 360.0*deg : ra_out[i];
		}
	}
}


void pixel_1d_to_pixel_2d_ref( int map_index,  int nx, int ny, int & out_x, 
    int & out_y)
{
	if (map_index >= nx * ny || map_index < 0){
		out_y = -1;
		out_x = -1;
	} else {
		out_x = map_index % nx;
		out_y = map_index / nx;
	}
}


PYBINDINGS("coordinateutils")
{
	using namespace boost::python;
	def("pixel_2d_to_pixel_1d", pixel_2d_to_pixel_1d);
	def("pixel_1d_to_pixel_2d", pixel_1d_to_pixel_2d);
}
