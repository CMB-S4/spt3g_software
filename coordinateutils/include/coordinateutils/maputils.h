#ifndef _COORDINATEUTILS_MAPUTILS_H
#define _COORDINATEUTILS_MAPUTILS_H

#include <vector>

#include <G3SkyMap.h>

void get_ra_dec_map_cpp(G3SkyMapConstPtr m, G3SkyMapPtr ra, G3SkyMapPtr dec);

void reproj_map(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin=1);

void reproj_fullsky_healpix_map(std::vector<double> in_map,
    G3SkyMapPtr out_map, bool nest=false, int rebin=1);

void maputils_pybindings(void);

#endif //_COORDINATEUTILS_MAPUTILS_H
