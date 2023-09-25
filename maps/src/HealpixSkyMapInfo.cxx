#include <cmath>

#include <pybindings.h>
#include <serialization.h>

#include <G3Units.h>

#include <maps/HealpixSkyMapInfo.h>
#include "chealpix.h"

static const double twothird = 2.0 / 3.0;
static const double twopi = 2 * M_PI;

HealpixSkyMapInfo::HealpixSkyMapInfo(size_t nside, bool nested, bool shifted)
{
	initialize(nside, nested, shifted);
}

HealpixSkyMapInfo::HealpixSkyMapInfo()
{
	initialize(0, false, false);
}

HealpixSkyMapInfo::HealpixSkyMapInfo(const HealpixSkyMapInfo & other)
{
	initialize(other.nside_, other.nested_, other.shifted_);
}

void
HealpixSkyMapInfo::initialize(size_t nside_or_npix, bool nested, bool shifted,
   bool is_npix)
{
	if (is_npix)
		SetNPix(nside_or_npix);
	else
		SetNSide(nside_or_npix);
	SetNested(nested);
	SetShifted(shifted);
}

template <class A> void
HealpixSkyMapInfo::load(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("nside", nside_);
	ar & make_nvp("nested", nested_);
	ar & make_nvp("shifted", shifted_);

	initialize(nside_, nested_, shifted_);
}

template <class A> void
HealpixSkyMapInfo::save(A &ar, unsigned v) const
{
	using namespace cereal;

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("nside", nside_);
	ar & make_nvp("nested", nested_);
	ar & make_nvp("shifted", shifted_);
}

std::string
HealpixSkyMapInfo::Description() const
{
	std::ostringstream os;

	os.precision(1);

	os << "Nside-" << nside_ << ", "
	   << (nested_ ? "nested" : "ring-ordered") << ", "
	   << "center alpha=" << (shifted_ ? 0 : 180) << " deg";

	return os.str();
}

bool
HealpixSkyMapInfo::IsCompatible(const HealpixSkyMapInfo & other) const
{
	return ((nside_ == other.nside_) && (nested_ == other.nested_));
}

double
HealpixSkyMapInfo::res() const
{
	return sqrt(4.0 * M_PI / npix_) * G3Units::rad;
}

void
HealpixSkyMapInfo::SetNPix(size_t npix)
{
	ssize_t nside = npix2nside64(npix);
	if (nside < 0)
		log_fatal("Invalid npix %zu", npix);

	SetNSide(nside);
}

void
HealpixSkyMapInfo::SetNSide(size_t nside)
{
	nside_ = nside;
	nring_ = 4 * nside_;
	size_t npface = nside_ * nside_;
	npix_ = 12 * npface;
	ncap_ = 2 * (npface - nside_);
	double fact2 = 4. / npix_;
	double fact1 = 2 * nside_ * fact2;

	rings_.clear();
	rings_.resize(nring_);

	for (size_t iring = 0; iring < nring_; iring++) {
		HealpixRingInfo ring;
		size_t northring = (iring > (nring_ / 2)) ? (nring_ - iring) : iring;

		if (northring < nside_) {
			double tmp = northring * northring * fact2;
			double costheta = 1 - tmp;
			double sintheta = sqrt(tmp * (2 - tmp));
			ring.theta = atan2(sintheta, costheta);
			ring.z = 1.0 - northring * northring * fact2;
			ring.npix = 4 * northring;
			ring.shift = 0.5;
			ring.pix0 = 2 * northring * (northring - 1);
		} else {
			ring.z = (2 * nside_ - northring) * fact1;
			ring.theta = acos(ring.z);
			ring.npix = 4 * nside_;
			ring.shift = (((northring - nside_) & 1) == 0) ? 0.5 : 0.0;
			ring.pix0 = ncap_ + (northring - nside_) * ring.npix;
		}

		if (northring != iring) {
			ring.z *= -1;
			ring.theta = M_PI - ring.theta;
			ring.pix0 = npix_ - ring.pix0 - ring.npix;
		}

		ring.delta = (M_PI_2 - ring.theta) * G3Units::rad;
		ring.dphi = twopi / ring.npix;
		ring.dalpha = ring.dphi * G3Units::rad;

		rings_[iring] = ring;
	}
}

void
HealpixSkyMapInfo::SetNested(bool nested)
{
	nested_ = nested;
}

void
HealpixSkyMapInfo::SetShifted(bool shifted)
{
	shifted_ = shifted;
}

std::pair<size_t, size_t>
HealpixSkyMapInfo::PixelToRing(size_t pixel) const
{
	static const std::pair<size_t, size_t> bad = {(size_t) -1, (size_t)-1};

	if (pixel >= npix_)
		return bad;

	if (nested_) {
		ssize_t pix = pixel;
		nest2ring64(nside_, pix, &pix);
		pixel = pix;
	}

	size_t iring;
	if (pixel < ncap_) /* North Polar cap */
		iring = (size_t)(0.5 * (1 + sqrt(1.5 + 2 * pixel)));
	else if (pixel < (npix_ - ncap_)) /* Equatorial region */
		iring = (size_t)((pixel - ncap_) / nring_ + nside_);
	else /* South Polar cap */
		iring = nring_ - (size_t)(0.5 * (1 + sqrt(2 * (npix_ - pixel) - 0.5)));
	if (iring >= nring_)
		return bad;

	const HealpixRingInfo & ring = rings_[iring];

	size_t ringpix = pixel - ring.pix0;
	if (ringpix >= ring.npix)
		return bad;

	if (shifted_)
		ringpix = (ringpix + ring.npix / 2) % ring.npix;

	return {iring, ringpix};
}

size_t
HealpixSkyMapInfo::RingToPixel(size_t iring, size_t ringpix) const
{
	if (iring >= nring_)
		return -1;
	const HealpixRingInfo & ring = rings_[iring];

	if (shifted_)
		ringpix = (ringpix + ring.npix / 2) % ring.npix;
	if (ringpix >= ring.npix)
		return -1;

	size_t pixel = ring.pix0 + ringpix;
	if (pixel >= npix_)
		return -1;

	if (nested_) {
		ssize_t pix = pixel;
		ring2nest64(nside_, pix, &pix);
		pixel = pix;
	}

	return pixel;
}

size_t
HealpixSkyMapInfo::AngleToPixel(double alpha, double delta) const
{
	if (std::isnan(alpha) || std::isnan(delta))
		return -1;

	double theta = (90 * G3Units::deg - delta) / G3Units::rad;
	if (theta < 0 || theta > M_PI)
		return -1;

	alpha /= G3Units::rad;
	while (alpha < 0)
		alpha += twopi;

	ssize_t outpix;
	if (nested_)
		ang2pix_nest64(nside_, theta, alpha, &outpix);
	else
		ang2pix_ring64(nside_, theta, alpha, &outpix);

	if (outpix < 0 || outpix >= npix_)
		return -1;

	return outpix;
}

std::vector<double>
HealpixSkyMapInfo::PixelToAngle(size_t pixel) const
{
	if (pixel < 0 || pixel >= npix_)
		return {0., 0.};

	double alpha, delta;
	if (nested_)
		pix2ang_nest64(nside_, pixel, &delta, &alpha);
	else
		pix2ang_ring64(nside_, pixel, &delta, &alpha);

	while(alpha > M_PI)
		alpha -= twopi;

	if (delta < 0 || delta > 180 * G3Units::deg)
		return {0., 0.};

	alpha *= G3Units::rad;
	delta = 90 * G3Units::deg - delta * G3Units::rad;

	return {alpha, delta};
}

size_t
HealpixSkyMapInfo::RebinPixel(size_t pixel, size_t scale) const
{
	ssize_t pix = pixel;

	if (!nested_)
		ring2nest64(nside_, pix, &pix);
	pix /= scale * scale;
	if (!nested_)
		nest2ring64(nside_ / scale, pix, &pix);

	return pix;
}

void
HealpixSkyMapInfo::GetRebinAngles(size_t pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas) const
{
	if (nside_ % scale != 0)
		log_fatal("Nside must be a multiple of rebinning scale");

	if (!nested_) {
		ssize_t pix = pixel;
		ring2nest64(nside_, pix, &pix);
		pixel = pix;
	}

	alphas = std::vector<double>(scale * scale);
	deltas = std::vector<double>(scale * scale);

	size_t nside_rebin = nside_ * scale;
	size_t pixmin = pixel * scale * scale;

	for (size_t i = 0; i < (scale * scale); i++) {
		size_t p = pixmin + i;
		double theta, phi;
		pix2ang_nest64(nside_rebin, p, &theta, &phi);

		while (phi > M_PI)
			phi -= twopi;

		alphas[i] = phi * G3Units::rad;
		deltas[i] = 90 * G3Units::deg - theta * G3Units::rad;
	}
}

size_t
HealpixSkyMapInfo::RingAbove(double z) const
{
	double za = fabs(z);
	if (za <= twothird)
		return nside_ * (2 - 1.5 * z);
	size_t iring = nside_ * sqrt(3 * (1 - za));
	return (z > 0) ? iring : (nring_ - iring - 1);
}

void
HealpixSkyMapInfo::GetInterpPixelsWeights(double alpha, double delta,
    std::vector<size_t> & pixels, std::vector<double> & weights) const
{
	alpha /= G3Units::rad;
	delta /= G3Units::rad;

	pixels = std::vector<size_t>(4, (size_t) -1);
	weights = std::vector<double>(4, 0);

	double theta = M_PI_2 - delta;
	g3_assert(!(theta < 0 || theta > M_PI));
	double z = cos(theta);

	double phi = (alpha < 0) ? (alpha + twopi) : alpha;

	size_t ir1 = RingAbove(z);
	size_t ir2 = ir1 + 1;

	double z1, z2;

	if (ir1 > 0) {
		const HealpixRingInfo & ring = rings_[ir1];
		z1 = ring.z;
		double tmp = phi / ring.dphi - ring.shift;
		ssize_t i1 = (tmp < 0) ? tmp - 1 : tmp;
		double w1 = (phi - (i1 + ring.shift) * ring.dphi) / ring.dphi;
		if (i1 < 0)
			i1 += ring.npix;
		ssize_t i2 = i1 + 1;
		if (i2 >= ring.npix)
			i2 -= ring.npix;
		pixels[0] = ring.pix0 + i1;
		pixels[1] = ring.pix0 + i2;
		weights[0] = 1 - w1;
		weights[1] = w1;
	}

	if (ir2 < nring_) {
		const HealpixRingInfo & ring = rings_[ir2];
		z2 = ring.z;
		double tmp = phi / ring.dphi - ring.shift;
		ssize_t i1 = (tmp < 0) ? tmp - 1 : tmp;
		double w1 = (phi - (i1 + ring.shift) * ring.dphi) / ring.dphi;
		if (i1 < 0)
			i1 += ring.npix;
		ssize_t i2 = i1 + 1;
		if (i2 >= ring.npix)
			i2 -= ring.npix;
		pixels[2] = ring.pix0 + i1;
		pixels[3] = ring.pix0 + i2;
		weights[2] = 1 - w1;
		weights[3] = w1;
	}

	if (ir1 == 0) {
		double wz = (z - 1.0) / (z2 - 1.0);
		weights[2] *= wz;
		weights[3] *= wz;
		double fac = (1 - wz) * 0.25;
		weights[0] = fac;
		weights[1] = fac;
		weights[2] += fac;
		weights[3] += fac;
		pixels[0] = (pixels[2] + 2) & 3;
		pixels[1] = (pixels[3] + 2) & 3;
	} else if (ir2 == nring_) {
		double wz = (z - z1) / (-1.0 - z1);
		weights[0] *= 1 - wz;
		weights[1] *= 1 - wz;
		double fac = wz * 0.25;
		weights[0] += fac;
		weights[1] += fac;
		weights[2] = fac;
		weights[3] = fac;
		pixels[2] = ((pixels[0] + 2) & 3) + npix_ - 4;
		pixels[3] = ((pixels[1] + 2) & 3) + npix_ - 4;
	} else {
		double wz = (z - z1) / (z2 - z1);
		weights[0] *= 1 - wz;
		weights[1] *= 1 - wz;
		weights[2] *= wz;
		weights[3] *= wz;
	}

	if (nested_) {
		for (size_t i = 0; i < pixels.size(); i++) {
			ssize_t pix = pixels[i];
			ring2nest64(nside_, pix, &pix);
			pixels[i] = pix;
		}
	}
}

std::vector<size_t>
HealpixSkyMapInfo::QueryDisc(double alpha, double delta, double radius) const
{
	auto pixels = std::vector<size_t>();

	radius /= G3Units::rad;
	if (radius >= M_PI) {
		for (ssize_t i = 0; i < npix_; i++)
			pixels.push_back(i);
		return pixels;
	}

	alpha /= G3Units::rad;
	delta /= G3Units::rad;

	double theta = M_PI_2 - delta;
	g3_assert(!(theta < 0 || theta > M_PI));
	double cosrad = cos(radius);
	double z = cos(theta);
	double xa = 1.0 / sqrt((1.0 - z) * (1.0 + z));
	double phi = (alpha < 0) ? (alpha + twopi) : alpha;

	double rlat1 = theta - radius;
	double zmax = cos(rlat1);
	size_t irmin = RingAbove(zmax) + 1;

	if ((rlat1 <= 0) && (irmin > 1)) {
		// north pole in the disk
		const HealpixRingInfo & ring = rings_[irmin - 1];
		for (ssize_t i = 0; i < ring.pix0 + ring.npix; i++)
			pixels.push_back(i);
	}

	double rlat2 = theta + radius;
	double zmin = cos(rlat2);
	size_t irmax = RingAbove(zmin);

	for (size_t iring = irmin; iring <= irmax; iring++) {
		const HealpixRingInfo & ring = rings_[iring];

		double x = (cosrad - ring.z * z) * xa;
		double ysq = 1 - ring.z * ring.z - x * x;
		if (ysq <= 0)
			continue;

		double dphi = atan2(sqrt(ysq), x);

		// highest pixel number in the ring
		ssize_t ipix2 = ring.pix0 + ring.npix - 1;

		ssize_t ip_lo = (ssize_t)floor((phi - dphi) / ring.dphi - ring.shift) + 1;
		ssize_t ip_hi = (ssize_t)floor((phi + dphi) / ring.dphi - ring.shift);
		if (ip_lo > ip_hi)
			continue;

		if (ip_hi >= ring.npix) {
			ip_lo -= ring.npix;
			ip_hi -= ring.npix;
		}
		if (ip_lo < 0) {
			for (ssize_t i = ring.pix0; i < ring.pix0 + ip_hi + 1; i++)
				pixels.push_back(i);
			for (ssize_t i = ring.pix0 + ip_lo + ring.npix; i < ipix2 + 1; i++)
				pixels.push_back(i);
		} else {
			for (ssize_t i = ring.pix0 + ip_lo; i < ring.pix0 + ip_hi + 1; i++)
				pixels.push_back(i);
		}
	}

	if ((rlat2 >= M_PI) && (irmax + 1 < nring_)) {
		// south pole in the disk
		const HealpixRingInfo & ring = rings_[irmax + 1];
		for (ssize_t i = ring.pix0; i < npix_; i++)
			pixels.push_back(i);
	}

	if (nested_) {
		for (size_t i = 0; i < pixels.size(); i++) {
			ssize_t pix = pixels[i];
			ring2nest64(nside_, pix, &pix);
			pixels[i] = pix;
		}
		std::sort(pixels.begin(), pixels.end());
	}

	return pixels;
}

G3_SPLIT_SERIALIZABLE_CODE(HealpixSkyMapInfo);
