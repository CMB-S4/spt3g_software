#ifndef _COORDINATEUTILS_FLATSKYPROJECTION_H
#define _COORDINATEUTILS_FLATSKYPROJECTION_H

#include <vector>

#include <G3Frame.h>
#include <G3Module.h>
#include <G3Timestream.h>
#include <G3Map.h>
#include <G3Vector.h>

#include <G3SkyMap.h>

#include <coordinateutils/FlatSkyMap.h>

typedef double angle_t;

int pixel_2d_to_pixel_1d( int x, int y, int nx, int ny);
std::vector<int> pixel_1d_to_pixel_2d( int map_index,  int nx, int ny);

void angle_to_pixel_1d_vec( const std::vector<double> & ras, 
    const std::vector<double> & decs,
    double ra0, double dec0,
    int n_x, int n_y,
    double reso_arcmin, double x_reso_arcmin,
    MapProjection proj,
    std::vector<int> & output_index   );

void pixel_1d_to_angle(const std::vector<int> & pix_ind,
    double ra0, double dec0, int n_x, int n_y,
    double pixel_res, double x_res,
    MapProjection proj, bool wrap_ra, 
    std::vector<double> & ra_out, std::vector<double> & dec_out );


#endif //#ifndef _COORDINATEUTILS_FLATSKYPROJECTION_H
