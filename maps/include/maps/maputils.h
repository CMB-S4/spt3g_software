#ifndef _MAPS_MAPUTILS_H
#define _MAPS_MAPUTILS_H

#include <vector>

#include <pybindings.h>

#include <maps/G3SkyMap.h>
#include <maps/FlatSkyMap.h>

G3SkyMapPtr GetMaskMap(G3SkyMapConstPtr m, bool zero_nans=false, bool zero_infs=false);

G3SkyMapPtr ApplyMask(G3SkyMapPtr m, G3SkyMapConstPtr mask, bool inverse=false);

void RemoveWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W, bool zero_nans = false);
void RemoveWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W,
    bool zero_nans = false);

void ApplyWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W);
void ApplyWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W);

boost::python::tuple GetRaDecMap(G3SkyMapConstPtr m);

G3SkyMapPtr GetRaDecMask(G3SkyMapConstPtr m, double ra_left, double ra_right,
    double dec_bottom, double dec_top);

void ReprojMap(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin=1, bool interp=false);

void FlattenPol(FlatSkyMapPtr Q, FlatSkyMapPtr U, double h=0.001, bool invert=false);

std::vector<double> GetMapStats(G3SkyMapConstPtr m, G3SkyMapConstPtr mask=NULL, int order=2,
    bool ignore_zeros=false, bool ignore_nans=false, bool ignore_infs=false);

double GetMapMedian(G3SkyMapConstPtr m, G3SkyMapConstPtr mask=NULL,
    bool ignore_zeros=false, bool ignore_nans=false, bool ignore_infs=false);

FlatSkyMapPtr ConvolveMap(FlatSkyMapConstPtr map, FlatSkyMapConstPtr kernel);

void maputils_pybindings(void);

#endif //_MAPS_MAPUTILS_H
