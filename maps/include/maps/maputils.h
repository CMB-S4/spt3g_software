#ifndef _MAPS_MAPUTILS_H
#define _MAPS_MAPUTILS_H

#include <vector>

#include <pybindings.h>

#include <maps/G3SkyMap.h>

boost::python::tuple get_ra_dec_map(G3SkyMapConstPtr m);

void reproj_map(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin=1, bool interp=false);

void maputils_pybindings(void);

#endif //_MAPS_MAPUTILS_H
