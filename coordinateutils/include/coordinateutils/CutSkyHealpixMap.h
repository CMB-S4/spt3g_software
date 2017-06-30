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
 *  sky pixel index
 *
 * If you apply the HealpixHitPix mapping of cut sky to full sky index and
 * then use the healpix routine pix_to_ang, it will return the appropriate
 * galactic coordinate system angles for where that pixel maps to. In galactic
 * coordinates, this allows a 1:1 mapping to Planck maps.
 *
 * The radec_2_pix routine accepts arguments in equatorial coordinates and then
 * performs an internal transform to galactic coordinates internally.
 *
 * Otherwise, the ang_2_pix and pix_2_ang methods maintain whatever coordinate
 * system is provided.
 */

class HealpixHitPix : public G3FrameObject {
public:
	HealpixHitPix(){};

	// Initialize with all the pixels in a flat sky map, potentially
	// transformed into a coordinate system given by coord_reference
	HealpixHitPix(const FlatSkyMap &flat_map, size_t nside, 
            bool is_nested, MapCoordReference coord_reference);

	// Initialize with a vector of pixel indices to include
	HealpixHitPix(const std::vector<uint64_t> &pixinds, size_t nside, 
	    bool is_nested, MapCoordReference coord_ref);

	long get_fullsky_index(long cutsky_index) const;
	long get_cutsky_index(long fullsky_index) const;
  
	bool is_nested() const { return is_nested_; }
	size_t get_nside() const { return nside_; }
	long get_ipixmax() const { return ipixmax_; }
	long get_ipixmin() const { return ipixmin_; }
	long get_total() const { return total_; }
	
	template <class A> void serialize(A &ar, unsigned u);
	std::string Description() const;

	MapCoordReference coord_ref;

public: // XXX: should be private
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
	// XXX: Document arguments
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

	void get_full_sized_healpix_map(G3VectorDouble & full_sized_map);
	
	void get_pointing_pixel(const std::vector<double> &alpha,
	    const std::vector<double> &delta, std::vector<int> &out_inds) const;
	
	/*
	 * The alpha and delta provided is in FK5 coordinates (equatorial)
	 * If the map is galactic translates the FK5 coordinates into 
	 * galactic and returns that pixel.
	 */
	long radec_2_pix(double alpha, double delta) const;

	// As above, but without the possible coordinate transform
	long ang_2_pix(double alpha, double delta) const;

	// Always returns RA, Dec
	void pix_2_ang(long pix, double & alpha, double & delta) const;

	// XXX: Document these interpolation functions
	void get_interpolated_weights(double alpha, double delta,
	    long pix[4], double weight[4]) const;
	double get_interp_precalc(long pix[4], double weight[4]) const;
	double get_interpolated_value(double alpha, double delta) const;
	
	template <class A> void serialize(A &ar, unsigned u);

	std::string Description() const;

	bool is_nested() const { return is_nested_; }
	size_t get_nside() const { return nside_; }

	HealpixHitPixPtr hitpix_; // XXX: should be renamed without _

private:
	SET_LOGGER("CutSkyHealpixMap");	
	
	size_t nside_;
	int is_nested_;
};

G3_POINTERS(CutSkyHealpixMap);
G3_SERIALIZABLE(CutSkyHealpixMap, 1);

#endif //_COORDINATEUTILS_CUTSKYHEALPIX_H

