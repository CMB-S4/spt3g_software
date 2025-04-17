#ifndef _MAPS_HEALPIXSKYMAPINFO_H
#define _MAPS_HEALPIXSKYMAPINFO_H

#include <vector>
#include <string>

#include <G3Frame.h>
#include <G3Logging.h>
#include <G3Quat.h>

class HealpixSkyMapInfo : public G3FrameObject {
public:
	HealpixSkyMapInfo(size_t nside_or_npix, bool nested=false, bool shifted=false,
	    bool is_npix=false);

	HealpixSkyMapInfo();
	HealpixSkyMapInfo(const HealpixSkyMapInfo& info);

	void initialize(size_t nside_or_npix, bool nested=false,
	    bool shifted=false, bool is_npix=false);

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;
	bool IsCompatible(const HealpixSkyMapInfo & other) const;
	std::string Description() const override;

	void SetNSide(size_t nside);
	void SetNPix(size_t npix);
	void SetNested(bool nested);
	void SetShifted(bool shifted);

	size_t nside() const { return nside_; }
	size_t npix() const { return npix_; }
	size_t nring() const { return nring_; }
	bool nested() const { return nested_; }
	bool shifted() const { return shifted_; }
	double res() const;

	std::pair<size_t, size_t> PixelToRing(size_t pix) const;
	size_t RingToPixel(size_t iring, size_t ringpix) const;

	std::vector<double> PixelToAngle(size_t pixel) const;
	size_t AngleToPixel(double alpha, double delta) const;
	Quat PixelToQuat(size_t pixel) const;
	size_t QuatToPixel(const Quat &q) const;

	size_t RebinPixel(size_t pixel, size_t scale) const;

	G3VectorQuat GetRebinQuats(size_t pixel, size_t scale) const;

	void GetInterpPixelsWeights(const Quat &q, std::vector<uint64_t> & pixels,
	    std::vector<double> & weights) const;

	std::vector<uint64_t> QueryDisc(const Quat &q, double radius) const;

private:
	// scheme
	size_t nside_;
	bool nested_;
	bool shifted_;

	// derived values
	size_t nring_, npix_, ncap_;

	struct HealpixRingInfo {
		size_t pix0, npix;
		double theta;
		double delta;
		double z;
		double shift;
		double dphi;
		double dalpha;
	};

	std::vector<HealpixRingInfo> rings_;

	size_t RingAbove(double z) const;

	SET_LOGGER("HealpixSkyMapInfo");
};

G3_POINTERS(HealpixSkyMapInfo);
G3_SPLIT_SERIALIZABLE(HealpixSkyMapInfo, 1);

#endif //#ifndef _MAPS_HEALPIXSKYMAPINFO_H
