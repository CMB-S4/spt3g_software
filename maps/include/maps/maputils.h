#ifndef _MAPS_MAPUTILS_H
#define _MAPS_MAPUTILS_H

#include <vector>

#include <pybindings.h>

#include <maps/G3SkyMap.h>
#include <maps/FlatSkyMap.h>

boost::python::tuple get_ra_dec_map(G3SkyMapConstPtr m);

void reproj_map(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin=1, bool interp=false);

void flatten_pol(FlatSkyMapPtr Q, FlatSkyMapPtr U, double h=0.001, bool invert=false);

void maputils_pybindings(void);

#endif //_MAPS_MAPUTILS_H
