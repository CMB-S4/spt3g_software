#ifndef _MAPS_MAPUTILS_H
#define _MAPS_MAPUTILS_H

#include <vector>

#include <pybindings.h>

#include <maps/G3SkyMap.h>
#include <maps/FlatSkyMap.h>

// Return a map that is set to one where the input map is non-zero.
// Optionally zero out nans and/or infs if zero_nans or zero_infs is true.
G3SkyMapPtr GetMaskMap(G3SkyMapConstPtr m, bool zero_nans=false, bool zero_infs=false);

// Zero out any values in the input map where the mask is zero
// If inverse is true, zero the input map where the mask is non-zero
void ApplyMask(G3SkyMapPtr m, G3SkyMapConstPtr mask, bool inverse=false);

// Divide out map weights from the sky maps
void RemoveWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W, bool zero_nans = false);
void RemoveWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W,
    bool zero_nans = false);

// Apply map weights to sky maps
void ApplyWeightsT(G3SkyMapPtr T, G3SkyMapWeightsConstPtr W);
void ApplyWeights(G3SkyMapPtr T, G3SkyMapPtr Q, G3SkyMapPtr U, G3SkyMapWeightsConstPtr W);

// Return maps of RA and Dec coordinates for each pixel in the input map
boost::python::tuple GetRaDecMap(G3SkyMapConstPtr m);

// Return a mask that is one for all pixels in the input map that are within the
// given ra and dec rectangular bounds, and zero otherwise
G3SkyMapPtr GetRaDecMask(G3SkyMapConstPtr m, double ra_left, double ra_right,
    double dec_bottom, double dec_top);

// Reproject the input map onto the output map grid, optionally
// oversampling the output pixels by the rebin factor, and/or sampling the input map
// with interpolation
void ReprojMap(G3SkyMapConstPtr in_map, G3SkyMapPtr out_map, int rebin=1, bool interp=false);

// Flatten or unflatten Q and U polarization maps using the polarization gradient across the map.
// The h parameter controls the pixel width over which a gradient is computed.
void FlattenPol(FlatSkyMapPtr Q, FlatSkyMapPtr U, double h=0.001, bool invert=false);

// Compute map moments up to fourth order (mean, var, skewness, kurtosis) of the input map,
// optionally excluding any pixels that are zero in the input mask, or zero/nan/inf in the input map
std::vector<double> GetMapStats(G3SkyMapConstPtr m, G3SkyMapConstPtr mask=NULL, int order=2,
    bool ignore_zeros=false, bool ignore_nans=false, bool ignore_infs=false);

// Compute the median of the input map, optionally excluding any pixels that are
// zero in the input mask, or zero/nan/inf in the input map
double GetMapMedian(G3SkyMapConstPtr m, G3SkyMapConstPtr mask=NULL,
    bool ignore_zeros=false, bool ignore_nans=false, bool ignore_infs=false);

// Compute the min and max of the values in the input map, optionally excluding
// any pixels that are zero in the input mask, or zero/nan/inf in the input map
std::vector<double> GetMapMinMax(G3SkyMapConstPtr m, G3SkyMapConstPtr mask=NULL,
    bool ignore_zeros=false, bool ignore_nans=false, bool ignore_infs=false);

// Compute the histogram of the input map pixels, grouping the values into bins
// defined by the array of bin edges, and optionally excluding any pixels that are
// zero in the input mask, or zero/nan/inf in the input map
std::vector<double> GetMapHist(G3SkyMapConstPtr m, const std::vector<double> &bin_edges,
    G3SkyMapConstPtr mask, bool ignore_zeros, bool ignore_nans, bool ignore_infs);

// Convolve the input flat sky map with a filter kernel
FlatSkyMapPtr ConvolveMap(FlatSkyMapConstPtr map, FlatSkyMapConstPtr kernel);

// Python bindings
void maputils_pybindings(void);

#endif //_MAPS_MAPUTILS_H
