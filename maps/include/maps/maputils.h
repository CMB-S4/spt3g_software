#ifndef _MAPS_MAPUTILS_H
#define _MAPS_MAPUTILS_H

#include <vector>

#include <pybindings.h>

#include <maps/G3SkyMap.h>
#include <maps/FlatSkyMap.h>

G3SkyMapPtr GetMaskMap(G3SkyMapConstPtr m);

void RemoveWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W, bool zero_nans = false);
void RemoveWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W,
    bool zero_nans = false);

void ApplyWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W);
void ApplyWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W);

boost::python::tuple GetRaDecMap(G3SkyMapConstPtr m);

void ReprojMap(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin=1, bool interp=false);

void FlattenPol(FlatSkyMapPtr Q, FlatSkyMapPtr U, double h=0.001, bool invert=false);

void maputils_pybindings(void);

#endif //_MAPS_MAPUTILS_H
