#ifndef _COORDINATEUTILS_CUTSKYHEALPIX_H
#define _COORDINATEUTILS_CUTSKYHEALPIX_H

#include <G3Frame.h>
#include <G3SkyMap.h>
#include <G3Module.h>

#include <cmath>
#include <cstdint>
#include <string>

#include <coordinateutils/FlatSkyMap.h>
#include <coordinateutils/chealpix.h>
#include <coordinateutils/wcs.h>
#include <coordinateutils/wcslib.h>

typedef map_info WCSMapInfo;
G3_POINTERS(WCSMapInfo);

/*
 * The Healpix library contains a mapping of full sky pixel index <-> angle on
 *  sky
 *
 * The HealpixHitPix class contains a mapping of full sky pixel index <-> cut
 *  sky pixel index.  It defines the area of sky for a cutsky healpix map
 *
 * CutSkyHealpixMap stores the actual map information.
 */


class CutSkyHealpixMap;

class HealpixHitPix : public G3FrameObject {
public:
	HealpixHitPix(){};
	
	HealpixHitPix(const FlatSkyMap &flat_map, size_t nside, 
	    bool is_nested, MapCoordReference coord_reference);

	// Initialize with a vector of pixel indices to include
	HealpixHitPix(const std::vector<uint64_t> &pixinds, size_t nside, 
	    bool is_nested, MapCoordReference coord_ref);

	long get_fullsky_index(long cutsky_index) const;
	long get_cutsky_index(long fullsky_index) const;
  
	bool is_nested() const { return is_nested_; }
	size_t get_nside() const { return nside_; }

	template <class A> void serialize(A &ar, unsigned u);
	std::string Description() const;

	MapCoordReference coord_ref;

private:
	friend class CutSkyHealpixMap;
	
	bool is_nested_;
	size_t nside_;	
	long ipixmin_;
	long ipixmax_;
	long total_; //total number of pixels in the map

	std::vector<uint64_t> pixinds_; // maps cutsky to full sky
	std::vector<uint64_t> full_to_cut_inds_; // maps full sky to cutsky, kind of
};

G3_POINTERS(HealpixHitPix);
G3_SERIALIZABLE(HealpixHitPix, 1);

class CutSkyHealpixMap : public G3SkyMap {
public:
	CutSkyHealpixMap(boost::python::object v, size_t full_sky_map_nside,
	    HealpixHitPixPtr hitpix, bool is_weighted = true,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    G3Timestream::TimestreamUnits u = G3Timestream::Kcmb,
	    int init_sym_group = 0);
	
	CutSkyHealpixMap(boost::python::object v, HealpixHitPixPtr hitpix, 
	    bool is_weighted = true,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    G3Timestream::TimestreamUnits u = G3Timestream::Kcmb);
	
	CutSkyHealpixMap(HealpixHitPixPtr hitpix, bool is_weighted = true,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    G3Timestream::TimestreamUnits u = G3Timestream::Kcmb);
	
	CutSkyHealpixMap();
	
	void get_full_sized_healpix_map(std::vector<double> & full_sized_map);
	
	void get_pointing_pixel(const std::vector<double> &alpha,
	    const std::vector<double> &delta, std::vector<int> &out_inds) const;
	
	std::vector<int> angles_to_pixels(const std::vector<double> & alphas, 
	    const std::vector<double> & deltas) const;
	
	std::vector<double> pixel_to_angle(size_t pix) const;
	
	void get_interpolated_weights(double alpha, double delta,
	    long pix[4], double weight[4]) const;
	double get_interp_precalc(long pix[4], double weight[4]) const;
	double get_interpolated_value(double alpha, double delta) const;
	
	template <class A> void serialize(A &ar, unsigned u);
	
	std::string Description() const;
	
	bool is_nested() const { return is_nested_; }
	size_t get_nside() const { return nside_; }
	HealpixHitPixPtr hitpix; 
	
private:
	long ang_2_pix_(double alpha, double delta) const;
	
	SET_LOGGER("CutSkyHealpixMap");	
	
	size_t nside_;
	int is_nested_;
};

G3_POINTERS(CutSkyHealpixMap);
G3_SERIALIZABLE(CutSkyHealpixMap, 1);

#endif //_COORDINATEUTILS_CUTSKYHEALPIX_H

