#ifndef _COORDINATEUTILS_CUTSKYHEALPIX_H
#define _COORDINATEUTILS_CUTSKYHEALPIX_H

#include <G3Frame.h>
#include <G3Module.h>

#include <cmath>
#include <cstdint>
#include <string>

#include <coordinateutils/G3SkyMap.h>
#include <coordinateutils/FlatSkyMap.h>
#include <coordinateutils/chealpix.h>

/*
 * The Healpix library contains a mapping of full sky pixel index <-> angle on
 *  sky
 *
 * The HealpixHitPix class contains a mapping of full sky pixel index <-> cut
 *  sky pixel index.  It defines the area of sky for a cutsky healpix map
 *
 * CutSkyHealpixMap stores the actual map information.
 */


typedef map_info HealpixMapInfo;
G3_POINTERS(HealpixMapInfo);

class CutSkyHealpixMap;

class HealpixHitPix : public G3FrameObject, public G3SkyMap {
public:
	HealpixHitPix(){};

        // Basic object, useful for manipulating fullsky pixels
	HealpixHitPix(size_t nside, bool is_nested,
	    MapCoordReference coord_ref);

	HealpixHitPix(const FlatSkyMap & in_map, size_t nside,
	    bool is_nested);
	HealpixHitPix(const CutSkyHealpixMap & in_map, size_t nside,
	    bool is_nested);

	// Initialize with a vector of pixel indices to include
	HealpixHitPix(const std::vector<int64_t> &pixinds, size_t nside,
	    bool is_nested, MapCoordReference coord_ref);

	long get_fullsky_index(long cutsky_index) const;
	long get_cutsky_index(long fullsky_index) const;

	bool is_nested() const { return is_nested_; }
	size_t get_nside() const { return nside_; }
	size_t get_size(void) const { return total_; }

	void pixels_from_map(const G3SkyMap & in_map);

	long angle_to_pixel(double alpha, double delta,
	    bool cutsky=true) const;
	std::vector<double> pixel_to_angle(long pixel,
	    bool cutsky=true) const;

        void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas,
	    bool cutsky=true) const;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights,
	    bool cutsky=true) const;

	template <class A> void serialize(A &ar, unsigned u);
	std::string Description() const;

	MapCoordReference coord_ref;

	std::vector<int64_t> pixinds_; // maps cutsky to full sky
private:
	friend class CutSkyHealpixMap;

	bool is_nested_;
	size_t nside_;
	long ipixmin_;
	long ipixmax_;
	long total_; //total number of pixels in the map (not counting the overflow)

	std::vector<int64_t> full_to_cut_inds_; // maps full sky to cutsky, kind of

	HealpixMapInfoPtr map_info_;
};

G3_POINTERS(HealpixHitPix);
G3_SERIALIZABLE(HealpixHitPix, 1);

class CutSkyHealpixMap : public G3FrameObject, public G3SkyMap {
public:
	CutSkyHealpixMap(boost::python::object v, size_t full_sky_map_nside,
	    HealpixHitPixPtr hitpix, bool is_weighted = true,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb,
	    int init_sym_group = 0);

	CutSkyHealpixMap(boost::python::object v, HealpixHitPixPtr hitpix,
	    bool is_weighted = true,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb);

	CutSkyHealpixMap(HealpixHitPixPtr hitpix, bool is_weighted = true,
	    G3SkyMap::MapPolType pol_type = G3SkyMap::None,
	    G3Timestream::TimestreamUnits u = G3Timestream::Tcmb);

	CutSkyHealpixMap();

	bool IsCompatible(const G3SkyMap & other) const override;

	std::vector<double> get_fullsky_map();

	size_t angle_to_pixel(double alpha, double delta) const override;
	std::vector<double> pixel_to_angle(size_t pixel) const override;

	void get_rebin_angles(long pixel, size_t scale,
	    std::vector<double> & alphas, std::vector<double> & deltas) const override;
	void get_interp_pixels_weights(double alpha, double delta,
	    std::vector<long> & pixels, std::vector<double> & weights) const override;

	G3SkyMapPtr rebin(size_t scale) const override;

	template <class A> void serialize(A &ar, unsigned u);
	virtual G3SkyMapPtr Clone(bool copy_data = true) const override;
	std::string Description() const override;

	bool is_nested() const { return is_nested_; }
	size_t get_nside() const { return nside_; }
	HealpixHitPixPtr hitpix;

private:
	SET_LOGGER("CutSkyHealpixMap");

	size_t nside_;
	int is_nested_;
};

G3_POINTERS(CutSkyHealpixMap);
G3_SERIALIZABLE(CutSkyHealpixMap, 1);

#endif //_COORDINATEUTILS_CUTSKYHEALPIX_H

