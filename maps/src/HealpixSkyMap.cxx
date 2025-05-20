#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <typeinfo>
#include <sys/types.h>

#include <G3Units.h>
#include <maps/HealpixSkyMap.h>
#include <maps/G3SkyMapMask.h>

#include "mapdata.h"

HealpixSkyMap::HealpixSkyMap(size_t nside, bool weighted, bool nested,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, bool shift_ra, G3SkyMap::MapPolConv pol_conv) :
      G3SkyMap(coord_ref, weighted, u, pol_type, pol_conv),
      info_(nside, nested, shift_ra),
      dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
}

HealpixSkyMap::HealpixSkyMap(const HealpixSkyMapInfo & info, bool weighted,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, G3SkyMap::MapPolConv pol_conv) :
      G3SkyMap(coord_ref, weighted, u, pol_type, pol_conv),
      info_(info),
      dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
}

HealpixSkyMap::HealpixSkyMap() :
    G3SkyMap(MapCoordReference::Local, false),
    info_(0, false, false),
    dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
}

HealpixSkyMap::HealpixSkyMap(const HealpixSkyMap & fm) :
    G3SkyMap(fm), info_(fm.info_), dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
	if (fm.dense_)
		dense_ = new std::vector<double>(*fm.dense_);
	else if (fm.ring_sparse_)
		ring_sparse_ = new SparseMapData<double>(*fm.ring_sparse_);
	else if (fm.indexed_sparse_)
		indexed_sparse_ = new std::unordered_map<uint64_t, double>(
		    *fm.indexed_sparse_);
}

HealpixSkyMap::~HealpixSkyMap()
{
	if (dense_)
		delete dense_;
	if (ring_sparse_)
		delete ring_sparse_;
	if (indexed_sparse_)
		delete indexed_sparse_;
}

template <class A> void
HealpixSkyMap::save(A &ar, unsigned v) const
{
	using namespace cereal;

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
	ar & make_nvp("info", info_);
	if (dense_) {
		ar & make_nvp("store", 3);
		ar & make_nvp("data", *dense_);
	} else if (ring_sparse_) {
		ar & make_nvp("store", 2);
		ar & make_nvp("data", *ring_sparse_);
	} else if (indexed_sparse_) {
		ar & make_nvp("store", 1);
		ar & make_nvp("data", *indexed_sparse_);
	} else {
		ar & make_nvp("store", 0);
	}
}

template <class A> void
HealpixSkyMap::load(A &ar, unsigned v)
{
	using namespace cereal;
	int store;
	uint32_t nside;
	bool nested, shifted;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));

	if (v < 3) {
		ar & make_nvp("nside", nside);
		ar & make_nvp("nested", nested);
	} else {
		ar & make_nvp("info", info_);
	}

	if (dense_) {
		delete dense_;
		dense_ = NULL;
	}
	if (ring_sparse_) {
		delete ring_sparse_;
		ring_sparse_ = NULL;
	}
	if (indexed_sparse_) {
		delete indexed_sparse_;
		indexed_sparse_ = NULL;
	}

	ar & make_nvp("store", store);
	switch (store) {
	case 3:
		dense_ = new std::vector<double>();
		ar & make_nvp("data", *dense_);
		break;
	case 2:
		ring_sparse_ = new SparseMapData<double>(1,1);
		ar & make_nvp("data", *ring_sparse_);
		break;
	case 1:
		indexed_sparse_ = new std::unordered_map<uint64_t, double>;
		ar & make_nvp("data", *indexed_sparse_);
		break;
	}

	if (v == 2)
		ar & make_nvp("shift_ra", shifted);
	else if (v == 1)
		shifted = false;

	if (v < 3)
		info_.initialize(nside, nested, shifted);
}

HealpixSkyMap::const_iterator::const_iterator(const HealpixSkyMap &map, bool begin) :
    map_(map)
{
	if (map_.dense_) {
		it_dense_ = begin ? map_.dense_->begin() : map_.dense_->end();
		index_ = begin ? 0 : map_.size();
	} else if (map_.ring_sparse_) {
		auto iter = begin ? map_.ring_sparse_->begin() : map_.ring_sparse_->end();
		j_ = iter.x;
		k_ = iter.y;
	} else if (map_.indexed_sparse_) {
		it_indexed_sparse_ = begin ? map_.indexed_sparse_->begin() : map_.indexed_sparse_->end();
	} else
		index_ = 0;

	set_value();
}

HealpixSkyMap::const_iterator::const_iterator(const HealpixSkyMap::const_iterator & iter) :
    index_(iter.index_), value_(iter.value_), map_(iter.map_)
{
	if (map_.dense_)
		it_dense_ = iter.it_dense_;
	else if (map_.ring_sparse_) {
		j_ = iter.j_;
		k_ = iter.k_;
	} else if (map_.indexed_sparse_)
		it_indexed_sparse_ = iter.it_indexed_sparse_;
}

void
HealpixSkyMap::const_iterator::set_value()
{
	if (map_.dense_) {
		value_.second = index_ < map_.size() ? *it_dense_ : 0;
	} else if (map_.ring_sparse_) {
		index_ = map_.info_.RingToPixel(j_, k_);
		if (index_ < 0 || index_ >= map_.size()) {
			index_ = map_.size();
			value_.second = 0;
		} else {
			value_.second = map_.ring_sparse_->at(j_, k_);
		}
	} else if (map_.indexed_sparse_) {
		if (it_indexed_sparse_ != map_.indexed_sparse_->end()) {
			index_ = it_indexed_sparse_->first;
			value_.second = it_indexed_sparse_->second;
		} else {
			index_ = map_.size();
			value_.second = 0;
		}
	}

	value_.first = index_;
}

HealpixSkyMap::const_iterator
HealpixSkyMap::const_iterator::operator++()
{
	if (map_.dense_) {
		++index_;
		++it_dense_;
	} else if (map_.ring_sparse_) {
		SparseMapData<double>::const_iterator iter(*map_.ring_sparse_, j_, k_);
		++iter;
		j_ = iter.x;
		k_ = iter.y;
	} else if (map_.indexed_sparse_)
		++it_indexed_sparse_;

	set_value();
	return *this;
}

void
HealpixSkyMap::ConvertToDense()
{
	if (dense_)
		return;

	dense_ = new std::vector<double>(size(), 0);

	if (ring_sparse_) {
		size_t idx;
		for (auto i = ring_sparse_->begin(); i != ring_sparse_->end();
		    i++) {
			idx = info_.RingToPixel(i.x, i.y);
			if (idx < 0 || idx >= size())
				continue;
			(*dense_)[idx] = (*i);
		}
		delete ring_sparse_;
		ring_sparse_ = NULL;
	} else if (indexed_sparse_) {
		for (auto i : *indexed_sparse_)
			(*dense_)[i.first] = i.second;
		delete indexed_sparse_;
		indexed_sparse_ = NULL;
	}
	/* otherwise leave it at zero since no data */
}

void
HealpixSkyMap::ConvertToRingSparse()
{
	if (ring_sparse_)
		return;

	ring_sparse_ = new SparseMapData<double>(info_.nring(), info_.nring());

	if (dense_) {
		auto *old_dense = dense_;
		dense_ = NULL;

		for (size_t i = 0; i < old_dense->size(); i++)
			if ((*old_dense)[i] != 0)
				(*this)[i] = (*old_dense)[i];
		
		delete old_dense;
	} else if (indexed_sparse_) {
		auto *oldis = indexed_sparse_;
		indexed_sparse_ = NULL;

		for (auto i : *oldis)
			if (i.second != 0)
				(*this)[i.first] = i.second;

		delete oldis;
	}
}

void
HealpixSkyMap::ConvertToIndexedSparse()
{
	if (indexed_sparse_)
		return;

	indexed_sparse_ = new std::unordered_map<uint64_t, double>;

	if (ring_sparse_) {
		size_t idx;
		for (auto i = ring_sparse_->begin(); i != ring_sparse_->end();
		    i++) {
			if ((*i) == 0)
				continue;
			idx = info_.RingToPixel(i.x, i.y);
			if (idx < 0 || idx >= size())
				continue;
			(*indexed_sparse_)[idx] = (*i);
		}
		delete ring_sparse_;
		ring_sparse_ = NULL;
	} else if (dense_) {
		for (size_t i = 0; i < dense_->size(); i++)
			if ((*dense_)[i] != 0)
				(*indexed_sparse_)[i] = (*dense_)[i];
		delete dense_;
		dense_ = NULL;
	}
}

void HealpixSkyMap::Compact(bool zero_nans)
{
	if (!dense_ && !indexed_sparse_ && !ring_sparse_)
		return;

	if (zero_nans) {
		for (auto i : *this) {
			if (i.second != i.second)
				(*this)[i.first] = 0;
		}
	}

	if (dense_) {
		ConvertToRingSparse();
	} else if (indexed_sparse_) {
		for (auto i = indexed_sparse_->begin(); i != indexed_sparse_->end(); ) {
			if (i->second == 0)
				i = indexed_sparse_->erase(i);
			else
				i++;
		}
	} else if (ring_sparse_) {
		ring_sparse_->compact();
	}
}

G3SkyMapPtr
HealpixSkyMap::Clone(bool copy_data) const
{
	if (copy_data)
		return std::make_shared<HealpixSkyMap>(*this);
	else
		return std::make_shared<HealpixSkyMap>(nside(), weighted,
		    nested(), coord_ref, units, pol_type, info_.shifted(), pol_conv);
}

double *
HealpixSkyMap::data()
{
	return dense_ == NULL ? nullptr : dense_->data();
}

double
HealpixSkyMap::at(size_t i) const
{
	if (i < 0 || i >= info_.npix())
		return 0;

	if (dense_)
		return (*dense_)[i];
	if (ring_sparse_) {
		auto ridx = info_.PixelToRing(i);
		return ring_sparse_->at(ridx.first, ridx.second);
	}
	if (indexed_sparse_) {
		try {
			return indexed_sparse_->at(i);
		} catch (const std::out_of_range& e) {
			return 0;
		}
	}

	return 0;
}

double &
HealpixSkyMap::operator [] (size_t i)
{
	g3_assert(!(i < 0 || i >= info_.npix()));

	if (dense_)
		return (*dense_)[i];

	if (indexed_sparse_)
		return (*indexed_sparse_)[i];

	if (!ring_sparse_)
		ring_sparse_ = new SparseMapData<double>(info_.nring(), info_.nring());

	auto ridx = info_.PixelToRing(i);
	return (*ring_sparse_)(ridx.first, ridx.second);
}


#define healpixskymap_arithmetic(op) \
G3SkyMap &HealpixSkyMap::operator op(const G3SkyMap &rhs) { \
	g3_assert(IsCompatible(rhs)); \
	g3_assert(units == rhs.units); \
	g3_assert(weighted == rhs.weighted); \
	try { \
		const HealpixSkyMap& b = dynamic_cast<const HealpixSkyMap &>(rhs); \
		if (dense_ || ring_sparse_ || indexed_sparse_) { \
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) { \
				for (auto i : b) { \
					if (i.second != 0) \
						(*this)[i.first] op i.second; \
				} \
		        } \
		} else { \
			if (b.dense_) { \
				ConvertToDense(); \
				for (size_t i = 0; i < dense_->size(); i++) \
					(*dense_)[i] op (*b.dense_)[i]; \
			} else if (b.ring_sparse_) { \
				SetShiftRa(b.info_.shifted()); \
				ConvertToRingSparse(); \
				(*ring_sparse_) op (*b.ring_sparse_); \
			} else if (b.indexed_sparse_) { \
				ConvertToIndexedSparse(); \
				for (auto i : *b.indexed_sparse_) \
					if (i.second != 0) \
						(*indexed_sparse_)[i.first] op i.second; \
			} \
		} \
		return *this; \
	} catch (const std::bad_cast& e) { \
		return G3SkyMap::operator op(rhs); \
	} \
}

healpixskymap_arithmetic(+=)
healpixskymap_arithmetic(-=)

G3SkyMap &HealpixSkyMap::operator *=(const G3SkyMap &rhs) {
	g3_assert(IsCompatible(rhs));
	if (units == G3Timestream::None)
		units = rhs.units;
	if (rhs.weighted and !(weighted))
		weighted = true;
	try {
		const HealpixSkyMap& b = dynamic_cast<const HealpixSkyMap &>(rhs);
		bool zero = false;
		if (dense_ || ring_sparse_ || indexed_sparse_) {
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (auto i : *this)
					(*this)[i.first] *= b.at(i.first);
			} else
				zero = true;
		} else {
			zero = true;
		}

		if (zero) {
			if (indexed_sparse_)
				delete indexed_sparse_;
			if (ring_sparse_)
				delete ring_sparse_;
			if (dense_)
				delete dense_;
			dense_ = NULL;
			ring_sparse_ = NULL;
			indexed_sparse_ = NULL;
		}

		return *this;
	} catch (const std::bad_cast& e) {
		return G3SkyMap::operator *=(rhs);
	}
	return *this;
}

G3SkyMap &
HealpixSkyMap::operator *=(const G3SkyMapMask &rhs)
{
	g3_assert(rhs.IsCompatible(*this));

	for (auto i : *this) {
		if (!rhs.at(i.first) && i.second != 0)
			(*this)[i.first] = 0;
	}

	return *this;
}

G3SkyMap &HealpixSkyMap::operator /=(const G3SkyMap &rhs) {
	g3_assert(IsCompatible(rhs));
	if (units == G3Timestream::None)
		units = rhs.units;
	if (rhs.weighted and !(weighted))
		weighted = true;
	try {
		const HealpixSkyMap& b = dynamic_cast<const HealpixSkyMap &>(rhs);
		bool zero = false;
		if (dense_) {
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (size_t i = 0; i < dense_->size(); i++)
					(*dense_)[i] /= b.at(i);
			} else
				zero = true;
		} else if (ring_sparse_) {
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (size_t i = 0; i < size(); i++) {
					double valb = b.at(i);
					double val = this->at(i);
					if (valb == 0 || valb != valb || val != 0)
						(*this)[i] /= valb;
				}
			} else
				zero = true;
		} else if (indexed_sparse_) {
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (size_t i = 0; i < size(); i++) {
					double val = this->at(i);
					double valb = b.at(i);
					if (valb == 0 || valb != valb || val != 0)
						(*indexed_sparse_)[i] /= valb;
				}
			} else
				zero = true;
		} else {
			if (b.dense_) {
				ConvertToDense();
				for (size_t i = 0; i < dense_->size(); i++)
					(*dense_)[i] /= (*b.dense_)[i];
			} else if (b.ring_sparse_) {
				SetShiftRa(b.info_.shifted());
				ConvertToRingSparse();
				for (size_t i = 0; i < size(); i++) {
					double valb = b.at(i);
					if (valb == 0 || valb != valb)
						(*this)[i] /= valb;
				}
			} else if (b.indexed_sparse_) {
				ConvertToIndexedSparse();
				for (size_t i = 0; i < size(); i++) {
					double valb = b.at(i);
					if (valb == 0 || valb != valb)
						(*indexed_sparse_)[i] /= valb;
				}
			} else
				zero = true;
		}

		if (zero) {
			ConvertToDense();
			for (size_t i = 0; i < dense_->size(); i++)
				(*dense_)[i] /= 0.0;
		}

		return *this;
	} catch (const std::bad_cast& e) {
		return G3SkyMap::operator /=(rhs);
	}
	return *this;
}

#define healpixskymap_addd(op) \
G3SkyMap & \
HealpixSkyMap::operator op(double b) \
{ \
	if (b == 0) \
		return *this; \
	if (!dense_) \
		ConvertToDense(); \
	for (size_t i = 0; i < dense_->size(); i++) \
		(*dense_)[i] op b; \
	return *this; \
}

healpixskymap_addd(+=)
healpixskymap_addd(-=)

G3SkyMap &
HealpixSkyMap::operator *=(double b)
{
	if (b == 0) {
		if (ring_sparse_)
			delete ring_sparse_;
		if (indexed_sparse_)
			delete indexed_sparse_;
		if (dense_)
			delete dense_;
		ring_sparse_ = NULL;
		indexed_sparse_ = NULL;
		dense_ = NULL;
		return *this;
	}

	if (dense_)
		for (size_t i = 0; i < dense_->size(); i++)
			(*dense_)[i] *= b;
	else if (ring_sparse_)
		(*ring_sparse_) *= b;
	else if (indexed_sparse_) {
		for (auto i : *indexed_sparse_)
			(*indexed_sparse_)[i.first] *= b;
	}

	return *this;
}

G3SkyMap &
HealpixSkyMap::operator /=(double b)
{
	if (b == 0)
		ConvertToDense();

	if (dense_)
		for (size_t i = 0; i < dense_->size(); i++)
			(*dense_)[i] /= b;
	else if (ring_sparse_)
		(*ring_sparse_) /= b;
	else if (indexed_sparse_) {
		for (auto i : *indexed_sparse_)
			(*indexed_sparse_)[i.first] /= b;
	}

	return *this;
}

std::string
HealpixSkyMap::Description() const
{
	std::ostringstream os;

	os.precision(1);

	os << info_.Description() << " in ";

	switch (coord_ref) {
	case Local:
		os << "local";
		break;
	case Equatorial:
		os << "equatorial";
		break;
	case Galactic:
		os << "galactic";
		break;
	default:
		os << "unknown";
	}
	switch (pol_conv) {
	case IAU:
		os << " IAU";
		break;
	case COSMO:
		os << " COSMO";
		break;
	default:
		break;
	}
	os << " coordinates (";

	switch (units) {
	case G3Timestream::Counts:
		os << "Counts";
		break;
	case G3Timestream::Current:
		os << "Current";
		break;
	case G3Timestream::Power:
		os << "Power";
		break;
	case G3Timestream::Resistance:
		os << "Resistance";
		break;
	case G3Timestream::Tcmb:
		os << "Tcmb";
		break;
	case G3Timestream::Angle:
		os << "Angle";
		break;
	case G3Timestream::Distance:
		os << "Distance";
		break;
	case G3Timestream::Voltage:
		os << "Voltage";
		break;
	case G3Timestream::Pressure:
		os << "Pressure";
		break;
	case G3Timestream::FluxDensity:
		os << "FluxDensity";
		break;
	case G3Timestream::Trj:
		os << "Trj";
		break;
	case G3Timestream::Frequency:
		os << "Frequency";
		break;
	default:
		break;
	}

	os << ", " << (weighted ? "" : "not ") << "weighted)";

	return os.str();
}

bool
HealpixSkyMap::IsCompatible(const G3SkyMap & other) const
{
	try {
		const HealpixSkyMap &healp = dynamic_cast<const HealpixSkyMap&>(other);
		return (G3SkyMap::IsCompatible(other) && info_.IsCompatible(healp.info_));
	} catch(const std::bad_cast& e) {
		return false;
	}
}

void
HealpixSkyMap::NonZeroPixels(std::vector<uint64_t> &indices,
    std::vector<double> &data) const
{
	indices.clear();
	data.clear();

	size_t npix = NpixAllocated();
	if (npix == 0)
		return;

	indices.reserve(npix);
	data.reserve(npix);

	for (auto i : *this) {
		if (i.second != 0) {
			indices.push_back(i.first);
			data.push_back(i.second);
		}
	}
}

void
HealpixSkyMap::ApplyMask(const G3SkyMapMask &mask, bool inverse)
{
	g3_assert(mask.IsCompatible(*this));

	for (auto i: *this)
		if (i.second != 0 && mask.at(i.first) == inverse)
			(*this)[i.first] = 0;
}


std::vector<size_t>
HealpixSkyMap::shape() const
{
	return {info_.npix()};
}

double
HealpixSkyMap::res() const
{
	return info_.res();
}

size_t HealpixSkyMap::NpixAllocated() const
{
	if (dense_)
		return dense_->size();
	if (ring_sparse_)
		return ring_sparse_->allocated();
	if (indexed_sparse_)
		return indexed_sparse_->size();
	return 0;
}

size_t HealpixSkyMap::NpixNonZero() const
{
	size_t sz = 0;
	if (dense_) {
		for (auto v: *dense_)
			if (v != 0)
				sz++;
	} else if (ring_sparse_) {
		return ring_sparse_->nonzero();
	} else if (indexed_sparse_) {
		for (auto i : *indexed_sparse_)
			if (i.second != 0)
				sz++;
	}
	return sz;
}

size_t
HealpixSkyMap::AngleToPixel(double alpha, double delta) const
{
	return info_.AngleToPixel(alpha, delta);
}

std::vector<double>
HealpixSkyMap::PixelToAngle(size_t pixel) const
{
	return info_.PixelToAngle(pixel);
}

size_t
HealpixSkyMap::QuatToPixel(const Quat &q) const
{
	return info_.QuatToPixel(q);
}

Quat
HealpixSkyMap::PixelToQuat(size_t pixel) const
{
	return info_.PixelToQuat(pixel);
}

G3VectorQuat
HealpixSkyMap::GetRebinQuats(size_t pixel, size_t scale) const
{
	return info_.GetRebinQuats(pixel, scale);
}

void
HealpixSkyMap::GetInterpPixelsWeights(const Quat &q, std::vector<uint64_t> & pixels,
    std::vector<double> & weights) const
{
	info_.GetInterpPixelsWeights(q, pixels, weights);
}

std::vector<uint64_t>
HealpixSkyMap::QueryDisc(const Quat &q, double radius) const
{
	return info_.QueryDisc(q, radius);
}

G3SkyMapPtr HealpixSkyMap::Rebin(size_t scale, bool norm) const
{
	if (nside() % scale != 0)
		log_fatal("Map nside must be a multiple of rebinning scale");

	if (scale <= 1)
		return Clone(true);

	HealpixSkyMapPtr out(new HealpixSkyMap(nside()/scale, weighted,
	    nested(), coord_ref, units, pol_type, info_.shifted(), pol_conv));

	if (dense_)
		out->ConvertToDense();
	else if (ring_sparse_)
		out->ConvertToRingSparse();
	else if (indexed_sparse_)
		out->ConvertToIndexedSparse();
	else
		return out;

	const double sqscal = norm ? scale * scale : 1.0;

	for (auto i : *this) {
		if (i.second == 0)
			continue;
		size_t ip = info_.RebinPixel(i.first, scale);
		(*out)[ip] += i.second / sqscal;
	}

	return out;
}

void HealpixSkyMap::SetShiftRa(bool shift)
{
	if (shift == info_.shifted())
		return;

	if (!IsRingSparse()) {
		info_.SetShifted(shift);
		return;
	}

	HealpixSkyMapInfo info(info_);
	info.SetShifted(shift);

	SparseMapData<double> *ring_sparse = new SparseMapData<double>(info_.nring(), info_.nring());
	for (auto i : *this) {
		if (i.second != 0) {
			auto ridx = info.PixelToRing(i.first);
			(*ring_sparse)(ridx.first, ridx.second) = i.second;
		}
	}

	delete ring_sparse_;
	info_.SetShifted(shift);
	ring_sparse_ = ring_sparse;
}

static void
HealpixSkyMap_fill(HealpixSkyMap &skymap, const py::cbuffer &v)
{
	auto info = v.request_contiguous();

	// Fall back to just 1-D
	if (info.ndim != 1)
		throw py::value_error("Only 1-D maps supported");

	size_t npix = info.shape[0];
	if (npix != skymap.size())
		log_fatal("Got array of shape (%zu,), expected (%zu,)", npix, skymap.size());

	skymap.ConvertToDense();
	double *d = &skymap[0];

	auto format = check_buffer_format(info.format);

	if (format == "d") {
		memcpy(d, info.ptr, skymap.size() * info.itemsize);
	} else if (format == "f") {
		for (size_t i = 0; i < skymap.size(); i++)
			d[i] = ((float *)info.ptr)[i];
	} else if (format == "i") {
		for (size_t i = 0; i < skymap.size(); i++)
			d[i] = ((int *)info.ptr)[i];
	} else if (format == "I") {
		for (size_t i = 0; i < skymap.size(); i++)
			d[i] = ((unsigned int *)info.ptr)[i];
	} else if (format == "l") {
		for (size_t i = 0; i < skymap.size(); i++)
			d[i] = ((long *)info.ptr)[i];
	} else if (format == "L") {
		for (size_t i = 0; i < skymap.size(); i++)
			d[i] = ((unsigned long *)info.ptr)[i];
	} else {
		throw py::type_error(std::string("Unknown type code ") + info.format);
	}
}

static HealpixSkyMapPtr
HealpixSkyMap_from_numpy(const py::array &v, bool weighted,
    bool nested, MapCoordReference coord_ref,
    G3Timestream::TimestreamUnits u, G3SkyMap::MapPolType pol_type,
    bool shift_ra, G3SkyMap::MapPolConv pol_conv)
{
	if (v.ndim() != 1)
			throw py::value_error("Only 1-D maps supported");

	HealpixSkyMapInfo info(v.shape(0), nested, shift_ra, true);
	HealpixSkyMapPtr skymap(new HealpixSkyMap(info, weighted,
	    coord_ref, u, pol_type, pol_conv));

	HealpixSkyMap_fill(*skymap, v);

	return skymap;
}

static HealpixSkyMapPtr
HealpixSkyMap_array_clone(const HealpixSkyMap &m, const py::array &v)
{
	auto skymap = std::dynamic_pointer_cast<HealpixSkyMap>(m.Clone(false));
	HealpixSkyMap_fill(*skymap, v);
	return skymap;
}

static void
HealpixSkyMap_fill_sparse(HealpixSkyMap &skymap, const py::array_t<int64_t> &index,
    const py::array_t<double> &data)
{
	if (index.size() != data.size())
		log_fatal("Index and data must have matching shapes.");

	if (index.ndim() != 1 || data.ndim() != 1)
		log_fatal("Index and data be 1D.");

	auto rindex = index.unchecked<1>();
	auto rdata = data.unchecked<1>();

	double phi_min = 2 * M_PI;
	double phi_max = 0;
	double phi_min_shift = 2 * M_PI;
	double phi_max_shift = 0;
	for (size_t i = 0; i < (size_t) index.size(); i++) {
		int64_t pix = unwrap_index(rindex(i), skymap.size());

		double ang = skymap.PixelToAngle(pix)[0];
		if (ang < 0)
			ang += 2 * M_PI * G3Units::rad;
		ang = fmod(ang, 2 * M_PI * G3Units::rad);
		if (ang < phi_min)
			phi_min = ang;
		if (ang > phi_max)
			phi_max = ang;
		ang = fmod(ang + M_PI * G3Units::rad, 2 * M_PI * G3Units::rad);
		if (ang < phi_min_shift)
			phi_min_shift = ang;
		if (ang > phi_max_shift)
			phi_max_shift = ang;
	}
	bool shift_ra = (phi_max - phi_min) > (phi_max_shift - phi_min_shift);

	skymap.SetShiftRa(shift_ra);
	skymap.ConvertToRingSparse();

	for (size_t i = 0; i < (size_t) index.size(); i++)
		skymap[rindex(i)] = rdata(i);
}

static HealpixSkyMapPtr
HealpixSkyMap_from_numpy_sparse(const py::array_t<int64_t> &index,
    const py::array_t<double> &data, size_t nside, bool weighted,
    bool nested, MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, G3SkyMap::MapPolConv pol_conv)
{
	HealpixSkyMapPtr skymap(new HealpixSkyMap(nside, weighted,
	    nested, coord_ref, u, pol_type, false, pol_conv));

	HealpixSkyMap_fill_sparse(*skymap, index, data);

	return skymap;
}

static void
HealpixSkyMap_setshiftra(HealpixSkyMap &m, bool v)
{
	m.SetShiftRa(v);
}

static void
HealpixSkyMap_setdense(HealpixSkyMap &m, bool v)
{
	if (v)
		m.ConvertToDense();
	else {
		throw py::value_error("Cannot set dense to False. "
		    "Set ringsparse or indexedsparse to True to convert from dense.");
	}
}

static void
HealpixSkyMap_setringsparse(HealpixSkyMap &m, bool v)
{
	if (v)
		m.ConvertToRingSparse();
	else {
		throw py::value_error("Cannot set ringsparse to False. "
		    "Set indexedsparse or dense to True to convert from ringsparse.");
	}
}

static void
HealpixSkyMap_setindexedsparse(HealpixSkyMap &m, bool v)
{
	if (v)
		m.ConvertToIndexedSparse();
	else {
		throw py::value_error("Cannot set indexedsparse to False. "
		    "Set ringsparse or dense to True to convert from indexedsparse.");
	}
}

static py::tuple
HealpixSkyMap_nonzeropixels(const HealpixSkyMap &m)
{
	auto i = std::vector<uint64_t>(); // XXX pointers?
	auto d = std::vector<double>();

	m.NonZeroPixels(i, d);

	return py::make_tuple(i, d);
}

static double
skymap_getitem(const G3SkyMap &skymap, ssize_t i)
{
	return skymap.at(unwrap_index(i, skymap.size()));
}

static void
skymap_setitem(G3SkyMap &skymap, ssize_t i, double val)
{
	skymap[unwrap_index(i, skymap.size())] = val;
}

static std::vector<double>
HealpixSkyMap_getitem_masked(const HealpixSkyMap &skymap, const G3SkyMapMask &m)
{
	g3_assert(m.IsCompatible(skymap));
	std::vector<double> out;

	for (auto i : skymap) {
		if (m.at(i.first))
			out.push_back(i.second);
	}

	return out;
}

static void
HealpixSkyMap_setitem_masked_scalar(HealpixSkyMap &skymap, const G3SkyMapMask &m,
    double dval)
{
	g3_assert(m.IsCompatible(skymap));

	for (auto i : skymap) {
		if (m.at(i.first))
			skymap[i.first] = dval;
	}
}

static void
HealpixSkyMap_setitem_masked(HealpixSkyMap &skymap, const G3SkyMapMask &m,
    const std::vector<double> &val)
{
	g3_assert(m.IsCompatible(skymap));

	if (val.size() != m.sum())
		throw py::value_error("Item dimensions do not match masked area");

	size_t j = 0;
	for (auto i : skymap) {
		if (m.at(i.first))
			skymap[i.first] = val[j++];
	}
}

static auto
HealpixSkyMap_buffer_info(HealpixSkyMap &m)
{
	m.ConvertToDense();

	return py::buffer_info(m.data(), sizeof(double), "d", 1,
	    {m.shape()[0]}, {sizeof(double)});
}

static void
HealpixSkyMap_setslice_1d(HealpixSkyMap &skymap, const py::slice &coords, const py::buffer &val)
{
	size_t start(0), stop(0), step(0), len(0);
	if (!coords.compute(skymap.size(), &start, &stop, &step, &len))
		throw py::error_already_set();
	if (start != 0 || stop != skymap.size())
		throw py::index_error("1D slicing not supported");

	HealpixSkyMap_fill(skymap, val);
}


G3_SPLIT_SERIALIZABLE_CODE(HealpixSkyMap);

PYBINDINGS("maps", scope)
{
	register_frameobject<HealpixSkyMap, G3SkyMap>(scope, "HealpixSkyMap",
	  py::buffer_protocol(),
	  "HealpixSkyMap is a G3SkyMap with the extra meta information about the "
	  "particular Healpix pixelization used.  In practice it behaves "
	  "(mostly) like a 1d numpy array.  If you find that you need numpy "
	  "functionality from a HealpixSkyMap, e.g. for slicing across the array, "
	  "you can access a numpy representation of the map using `np.asarray(m)`. "
	  "This does not copy the data, so any changes to the resulting array will "
	  "affect the data stored in the map.")
	    .def(py::init<>())
	    .def(py::init<size_t, bool, bool, MapCoordReference,
	        G3Timestream::TimestreamUnits, G3SkyMap::MapPolType, bool,
	        G3SkyMap::MapPolConv>(),
	        py::arg("nside"),
	        py::arg("weighted") = true,
	        py::arg("nested") = false,
	        py::arg("coord_ref") = MapCoordReference::Equatorial,
	        py::arg("units") = G3Timestream::Tcmb,
	        py::arg("pol_type") = G3SkyMap::None,
	        py::arg("shift_ra") = false,
	        py::arg("pol_conv") = G3SkyMap::ConvNone,
	        "Instantiate a HealpixSkyMap with given nside")
	    .def(py::init(&HealpixSkyMap_from_numpy_sparse),
	        py::arg("index"),
	        py::arg("data"),
	        py::arg("nside"),
	        py::arg("weighted") = true,
	        py::arg("nested") = false,
	        py::arg("coord_ref") = MapCoordReference::Equatorial,
	        py::arg("units") = G3Timestream::Tcmb,
	        py::arg("pol_type") = G3SkyMap::None,
	        py::arg("pol_conv") = G3SkyMap::ConvNone,
	        "Instantiate a sparse HealpixSkyMap from existing index and data "
	        "arrays corresponding to the given nside")
	    .def(py::init(&HealpixSkyMap_from_numpy),
	        py::arg("data"),
	        py::arg("weighted") = true,
	        py::arg("nested") = false,
	        py::arg("coord_ref") = MapCoordReference::Equatorial,
	        py::arg("units") = G3Timestream::Tcmb,
	        py::arg("pol_type") = G3SkyMap::None,
	        py::arg("shift_ra") = false,
	        py::arg("pol_conv") = G3SkyMap::ConvNone,
	        "Instantiate a dense Healpix map from an existing dense map.")

	    .def("array_clone", &HealpixSkyMap_array_clone, py::arg("array"),
	       "Return a map of the same type, populated with a copy of the input "
	       "numpy array")
	    .def_buffer(&HealpixSkyMap_buffer_info)
	    .def_property_readonly("nside", &HealpixSkyMap::nside, "Healpix resolution parameter")
	    .def_property_readonly("res", &HealpixSkyMap::res, "Map resolution in angular units")
	    .def_property_readonly("nested", &HealpixSkyMap::nested,
		"True if pixel ordering is nested, False if ring-ordered")
	    .def_property("shift_ra", &HealpixSkyMap::IsRaShifted, HealpixSkyMap_setshiftra,
		"True if the ringsparse representation of the map is stored "
		"with the rings centered at ra = 0 deg, rather than ra = 180 deg.")
	    .def_property("dense", &HealpixSkyMap::IsDense, HealpixSkyMap_setdense,
		"True if the map is stored with all elements, False otherwise. "
		"If set to True, converts the map to a dense representation." )
	    .def_property("ringsparse", &HealpixSkyMap::IsRingSparse, HealpixSkyMap_setringsparse,
		"True if the map is stored as a dense 2D region using ring "
		"ordering (analogous to FlatSkyMap's sparse mode). "
		"Ring-sparsity is efficient for dense blocks on a ring-ordered "
		"map (e.g. a continuous sky region), but is inefficient "
		"otherwise (e.g. nested pixel ordering or discontinous coverage). "
		"If set to True, converts the map to this representation." )
	    .def_property("indexedsparse", &HealpixSkyMap::IsIndexedSparse,
		HealpixSkyMap_setindexedsparse,
		"True if the map is stored as a list of non-zero pixels "
		"and values. More efficient than ring-sparse for maps with "
		"holes or very small filling factors. "
		"If set to True, converts the map to this representation." )

	    .def("__getitem__", &skymap_getitem)
	    .def("__setitem__", &skymap_setitem)
	    .def("__setitem__", HealpixSkyMap_setslice_1d)
	    .def("__getitem__", HealpixSkyMap_getitem_masked)
	    .def("__setitem__", HealpixSkyMap_setitem_masked)
	    .def("__setitem__", HealpixSkyMap_setitem_masked_scalar)

	    .def("nonzero_pixels", &HealpixSkyMap_nonzeropixels,
		"Returns a list of the indices of the non-zero pixels in the "
		"map and a list of the values of those non-zero pixels.")
	;
}

