#include <pybindings.h>
#include <serialization.h>
#include <typeinfo>

#include <coordinateutils/HealpixSkyMap.h>
#include <coordinateutils/chealpix.h>

#include "mapdata.h"

HealpixSkyMap::HealpixSkyMap(size_t nside, bool is_weighted, bool is_nested,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, bool shift_ra) :
      G3SkyMap(coord_ref, is_weighted, u, pol_type),
      nside_(nside), is_nested_(is_nested), dense_(NULL), ring_sparse_(NULL),
      indexed_sparse_(NULL), shift_ra_(shift_ra)
{
	ring_info_ = init_map_info(nside, is_nested, shift_ra, 1);
	npix_ = nside2npix(nside);
}

HealpixSkyMap::HealpixSkyMap(boost::python::object v, bool is_weighted,
    bool is_nested, MapCoordReference coord_ref,
    G3Timestream::TimestreamUnits u, G3SkyMap::MapPolType pol_type) :
      G3SkyMap(coord_ref, is_weighted, u, pol_type), is_nested_(is_nested),
      dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL), shift_ra_(false)
{
	Py_buffer view;

	if (boost::python::extract<size_t>(v).check()) {
		// size_t from Python is also a bp::object,
		// so a Python caller intending to call the above
		// constructor can get here by accident since
		// the signatures are degenerate. Handle the
		// confusion gracefully.
		nside_ = boost::python::extract<size_t>(v)();
		ring_info_ = init_map_info(nside_, is_nested_, shift_ra_, 1);
		npix_ = nside2npix(nside_);
		return;
	}

	if (PyTuple_Check(v.ptr()) && PyTuple_Size(v.ptr()) == 3) {
		// One option is that we got passed a tuple of numpy
		// arrays: first indices, next data, next nside.
		Py_buffer indexview, dataview;
#if PY_MAJOR_VERSION < 3
		if (PyInt_Check(PyTuple_GetItem(v.ptr(), 2))) {
			nside_ = PyInt_AsSsize_t(PyTuple_GetItem(v.ptr(), 2));
		} else 
#endif
		if (PyLong_Check(PyTuple_GetItem(v.ptr(), 2))) {
#if PY_MAJOR_VERSION < 3
			nside_ = PyLong_AsUnsignedLong(PyTuple_GetItem(v.ptr(), 2));
#else
			nside_ = PyLong_AsSize_t(PyTuple_GetItem(v.ptr(), 2));
#endif
		} else {
			PyErr_SetString(PyExc_TypeError,
			    "Third tuple element for sparse maps needs to be "
			    "nside");
			throw bp::error_already_set();
		}

		if (PyObject_GetBuffer(PyTuple_GetItem(v.ptr(), 0), &indexview,
		    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1)
			throw bp::error_already_set();

		if (PyObject_GetBuffer(PyTuple_GetItem(v.ptr(), 1), &dataview,
		    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
			PyBuffer_Release(&indexview);
			throw bp::error_already_set();
		}

		if (indexview.ndim != 1 || dataview.ndim != 1) {
			PyBuffer_Release(&indexview);
			PyBuffer_Release(&dataview);
			log_fatal("Only 1-D maps supported");
		}

		if (indexview.len/indexview.itemsize !=
		    dataview.len/dataview.itemsize) {
			PyBuffer_Release(&indexview);
			PyBuffer_Release(&dataview);
			log_fatal("Index and data must have matching shapes.");
		}

		if (strcmp(indexview.format, "l") != 0 &&
		    strcmp(indexview.format, "L") != 0) {
			PyBuffer_Release(&indexview);
			PyBuffer_Release(&dataview);
			log_fatal("Indices must be (long) integers.");
		}

		if (strcmp(dataview.format, "d") != 0) {
			PyBuffer_Release(&indexview);
			PyBuffer_Release(&dataview);
			log_fatal("Data must be double-precision (float64).");
		}

		size_t sz = indexview.len / indexview.itemsize;
		double phi_min = 2 * M_PI;
		double phi_max = 0;
		double phi_min_shift = 2 * M_PI;
		double phi_max_shift = 0;
		for (size_t i = 0; i < sz; i++) {
			unsigned long pix = ((unsigned long *)indexview.buf)[i];
			double ang = pixel_to_angle(pix)[0];
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
		shift_ra_ = (phi_max - phi_min) > (phi_max_shift - phi_min_shift);

		ring_info_ = init_map_info(nside_, is_nested_, shift_ra_, 1);
		npix_ = nside2npix(nside_);
		ring_sparse_ = new SparseMapData(ring_info_->nring, ring_info_->nring);

		for (size_t i = 0; i < sz; i++)
			(*this)[((unsigned long *)indexview.buf)[i]]=
			    ((double *)dataview.buf)[i];
		PyBuffer_Release(&indexview);
		PyBuffer_Release(&dataview);
		return;
	}

	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		// Fall back to just 1-D
		if (view.ndim != 1)
			log_fatal("Only 1-D maps supported");

		nside_ = npix2nside(view.shape[0]);
		dense_ = new std::vector<double>(view.shape[0]);

		double *d = &(*dense_)[0];
		if (strcmp(view.format, "d") == 0) {
			memcpy(d, view.buf, view.len);
		} else if (strcmp(view.format, "f") == 0) {
			for (size_t i = 0; i < view.len/sizeof(float); i++)
				d[i] = ((float *)view.buf)[i];
		} else if (strcmp(view.format, "i") == 0) {
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				d[i] = ((int *)view.buf)[i];
		} else if (strcmp(view.format, "I") == 0) {
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				d[i] = ((unsigned int *)view.buf)[i];
		} else if (strcmp(view.format, "l") == 0) {
			for (size_t i = 0; i < view.len/sizeof(long); i++)
				d[i] = ((unsigned long *)view.buf)[i];
		} else {
			log_fatal("Unknown type code %s", view.format);
		}
		PyBuffer_Release(&view);

		ring_info_ = init_map_info(nside_, is_nested_, shift_ra_, 1);
		npix_ = nside2npix(nside_);

		return;
	}

	throw bp::error_already_set();
}

HealpixSkyMap::HealpixSkyMap() :
    G3SkyMap(MapCoordReference::Local, false), nside_(0), is_nested_(false),
    dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL), shift_ra_(false)
{
	ring_info_ = init_map_info(nside_, is_nested_, shift_ra_, 1);
	npix_ = nside2npix(nside_);
}

HealpixSkyMap::HealpixSkyMap(const HealpixSkyMap & fm) :
    G3SkyMap(fm), nside_(fm.nside_), is_nested_(fm.is_nested_),
    dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL), shift_ra_(fm.shift_ra_)
{
	if (fm.dense_)
		dense_ = new std::vector<double>(*fm.dense_);
	else if (fm.ring_sparse_)
		ring_sparse_ = new SparseMapData(*fm.ring_sparse_);
	else if (fm.indexed_sparse_)
		indexed_sparse_ = new std::unordered_map<uint64_t, double>(
		    *fm.indexed_sparse_);
	ring_info_ = init_map_info(nside_, is_nested_, shift_ra_, 1);
	npix_ = nside2npix(nside_);
}

HealpixSkyMap::~HealpixSkyMap()
{
	if (dense_)
		delete dense_;
	if (ring_sparse_)
		delete ring_sparse_;
	if (indexed_sparse_)
		delete indexed_sparse_;

	free_map_info(ring_info_);
}

template <class A> void
HealpixSkyMap::save(A &ar, unsigned v) const
{
	using namespace cereal;

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
	ar & make_nvp("nside", nside_);
	ar & make_nvp("nested", is_nested_);
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
	ar & make_nvp("shift_ra", shift_ra_);
}

template <class A> void
HealpixSkyMap::load(A &ar, unsigned v)
{
	using namespace cereal;
	int store;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
	ar & make_nvp("nside", nside_);
	ar & make_nvp("nested", is_nested_);

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
		ring_sparse_ = new SparseMapData(1,1);
		ar & make_nvp("data", *ring_sparse_);
		break;
	case 1:
		indexed_sparse_ = new std::unordered_map<uint64_t, double>;
		ar & make_nvp("data", *indexed_sparse_);
		break;
	}

	if (v > 1)
		ar & make_nvp("shift_ra", shift_ra_);
	else
		shift_ra_ = false;

	free_map_info(ring_info_);
	ring_info_ = init_map_info(nside_, is_nested_, shift_ra_, 1);
	npix_ = nside2npix(nside_);
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
		long idx;
		get_pixel_index(map_.ring_info_, j_, k_, &idx);
		index_ = idx;
		value_.second = map_.ring_sparse_->at(j_, k_);
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
		SparseMapData::const_iterator iter(*map_.ring_sparse_, j_, k_);
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

	dense_ = new std::vector<double>(npix_, 0);

	if (ring_sparse_) {
		long idx;
		for (auto i = ring_sparse_->begin(); i != ring_sparse_->end(); ++i) {
			get_pixel_index(ring_info_, i.x, i.y, &idx);
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
	if (ring_sparse_ || is_nested_)
		return;

	ring_sparse_ = new SparseMapData(ring_info_->nring, ring_info_->nring);

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
		long idx;
		for (auto i = ring_sparse_->begin(); i != ring_sparse_->end(); ++i) {
			get_pixel_index(ring_info_, i.x, i.y, &idx);
			if ((*i) != 0)
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

G3SkyMapPtr
HealpixSkyMap::Clone(bool copy_data) const
{
	if (copy_data)
		return boost::make_shared<HealpixSkyMap>(*this);
	else
		return boost::make_shared<HealpixSkyMap>(nside_, is_weighted,
		    is_nested_, coord_ref, units, pol_type);
}

double
HealpixSkyMap::operator [] (size_t i) const
{
	if (i < 0 || i >= npix_)
		return 0;

	if (dense_)
		return (*dense_)[i];
	if (ring_sparse_) {
		long j, k;
		if (get_ring_index(ring_info_, i, &j, &k))
			return 0;
		return ring_sparse_->at(j, k);
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
	assert(i >= 0);
	assert(i < npix_);

	if (dense_)
		return (*dense_)[i];

	if (indexed_sparse_)
		return (*indexed_sparse_)[i];

	if (!ring_sparse_)
		ring_sparse_ = new SparseMapData(ring_info_->nring, ring_info_->nring);

	long j, k;
	int check = get_ring_index(ring_info_, i, &j, &k);
	assert(!check);
	return (*ring_sparse_)(j, k);
}


#define healpixskymap_arithmetic(op) \
G3SkyMap &HealpixSkyMap::operator op(const G3SkyMap &rhs) { \
	assert(IsCompatible(rhs)); \
	try { \
		const HealpixSkyMap& b = dynamic_cast<const HealpixSkyMap &>(rhs); \
		if (dense_) { \
			if (b.dense_) { \
				for (size_t i = 0; i < dense_->size(); i++) \
					(*dense_)[i] op (*b.dense_)[i]; \
			} else if (b.ring_sparse_) { \
				for (auto i : b) \
					(*dense_)[i.first] op i.second; \
			} else if (b.indexed_sparse_) { \
				for (auto i : *b.indexed_sparse_) \
					(*dense_)[i.first] op i.second; \
			} \
		} else if (ring_sparse_) { \
			if (b.dense_) { \
				for (size_t i = 0; i < b.dense_->size(); i++) { \
					double val = (*b.dense_)[i]; \
					if (val != 0) \
						(*this)[i] op val; \
				} \
			} else if (b.ring_sparse_) { \
				(*ring_sparse_) op (*b.ring_sparse_); \
			} else if (b.indexed_sparse_) { \
				for (auto i : *b.indexed_sparse_) \
					(*this)[i.first] op i.second; \
			} \
		} else if (indexed_sparse_) { \
			if (b.dense_) { \
				for (size_t i = 0; i < b.dense_->size(); i++) { \
					double val = (*b.dense_)[i]; \
					if (val != 0) \
						(*indexed_sparse_)[i] op val; \
				} \
			} else if (b.ring_sparse_) { \
				for (auto i : b) \
					(*indexed_sparse_)[i.first] op i.second; \
			} else if (b.indexed_sparse_) { \
				for (auto i : *b.indexed_sparse_) \
					(*indexed_sparse_)[i.first] op i.second; \
			} \
		} else { \
			if (b.dense_) { \
				ConvertToDense(); \
				for (size_t i = 0; i < dense_->size(); i++) \
					(*dense_)[i] op (*b.dense_)[i]; \
			} else if (b.ring_sparse_) { \
				ConvertToRingSparse(); \
				(*ring_sparse_) op (*b.ring_sparse_); \
			} else if (b.indexed_sparse_) { \
				ConvertToIndexedSparse(); \
				for (auto i : *b.indexed_sparse_) \
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
	assert(IsCompatible(rhs));
	try {
		const HealpixSkyMap& b = dynamic_cast<const HealpixSkyMap &>(rhs);
		bool zero = false;
		if (dense_ || ring_sparse_ || indexed_sparse_) {
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (auto i : *this)
					(*this)[i.first] *= b[i.first];
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

G3SkyMap &HealpixSkyMap::operator /=(const G3SkyMap &rhs) {
	assert(IsCompatible(rhs));
	try {
		const HealpixSkyMap& b = dynamic_cast<const HealpixSkyMap &>(rhs);
		bool zero = false;
		if (dense_) {
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (size_t i = 0; i < dense_->size(); i++)
					(*dense_)[i] /= b[i];
			} else
				zero = true;
		} else if (ring_sparse_) {
			long i = 0;
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (long j = 0; j < ring_info_->nring; j++) {
					for (long k = 0; k < ring_info_->rings[j].ringpix; k++) {
						double val = ring_sparse_->at(j, k);
						double valb = b[i];
						if (valb == 0 || val != 0)
							(*ring_sparse_)(j, k) /= valb;
						i++;
					}
				}
			} else
				zero = true;
		} else if (indexed_sparse_) {
			if (b.dense_ || b.ring_sparse_ || b.indexed_sparse_) {
				for (size_t i = 0; i < size(); i++) {
					double val = (*this)[i];
					double valb = b[i];
					if (valb == 0 || val != 0)
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
				ConvertToRingSparse();
				for (long j = 0; j < ring_info_->nring; j++) {
					for (long k = 0; k < ring_info_->rings[j].ringpix; k++)
						(*ring_sparse_)(j, k) /= b.ring_sparse_->at(j, k);
				}
			} else if (b.indexed_sparse_) {
				ConvertToIndexedSparse();
				for (size_t i = 0; i < size(); i++)
					(*indexed_sparse_)[i] /= b[i];
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

	os << "Nside-" << nside_ << " map in ";

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
	os << " coordinates ";

	switch (units) {
	case G3Timestream::Counts:
		os << " (Counts)";
		break;
	case G3Timestream::Current:
		os << " (Current)";
		break;
	case G3Timestream::Power:
		os << " (Power)";
		break;
	case G3Timestream::Resistance:
		os << " (Resistance)";
		break;
	case G3Timestream::Tcmb:
		os << " (Tcmb)";
		break;
	default:
		break;
	}

	return os.str();
}

bool
HealpixSkyMap::IsCompatible(const G3SkyMap & other) const
{
	try {
		const HealpixSkyMap &healp = dynamic_cast<const HealpixSkyMap&>(other);
		return (G3SkyMap::IsCompatible(other) &&
			nside_ == healp.nside_);
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

	size_t npix = npix_allocated();
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


std::vector<size_t>
HealpixSkyMap::shape() const
{
	std::vector<size_t> shape(1);
	shape[0] = npix_;
	return shape;
}

size_t HealpixSkyMap::npix_allocated() const
{
	if (dense_)
		return dense_->size();
	if (ring_sparse_)
		return ring_sparse_->nonzero();
	if (indexed_sparse_)
		return indexed_sparse_->size();
	return 0;
}

size_t
HealpixSkyMap::angle_to_pixel(double alpha, double delta) const
{
	long outpix = -1;
	double theta = (90 * G3Units::deg - delta) / G3Units::rad;

	alpha /= G3Units::rad;

	if ( std::isnan(theta) || std::isnan(alpha) ) {
		return -1;
	}

	if (is_nested_)
		ang2pix_nest(nside_, theta, alpha, &outpix);
	else
		ang2pix_ring(nside_, theta, alpha, &outpix);

	return outpix;
}

std::vector<double>
HealpixSkyMap::pixel_to_angle(size_t pixel) const
{
	double alpha, delta;
	if (is_nested_)
		pix2ang_nest(nside_, pixel, &delta, &alpha);
	else
		pix2ang_ring(nside_, pixel, &delta, &alpha);

	alpha *= G3Units::rad;
	delta = 90 * G3Units::deg - delta * G3Units::rad;

	std::vector<double> out = {alpha, delta};
	return out;
}

void HealpixSkyMap::get_rebin_angles(long pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas) const
{
	if (nside_ % scale != 0)
		log_fatal("Nside must be a multiple of rebinning scale");

	if (!is_nested_)
		ring2nest(nside_, pixel, &pixel);

	alphas = std::vector<double>(scale * scale);
	deltas = std::vector<double>(scale * scale);

	size_t nside_rebin = nside_ * scale;
	long pixmin = pixel * scale * scale;
	for (size_t i = 0; i < (scale * scale); i++) {
		long p = pixmin + i;
		double theta, phi;
		pix2ang_nest(nside_rebin, p, &theta, &phi);
		alphas[i] = phi * G3Units::rad;
		deltas[i] = 90 * G3Units::deg - theta * G3Units::rad;
	}
}

void
HealpixSkyMap::get_interp_pixels_weights(double alpha, double delta,
    std::vector<long> & pixels, std::vector<double> & weights) const
{
	alpha /= G3Units::rad;
	delta /= G3Units::rad;

	double theta = M_PI/2.0 - delta;
	double w[4];
	long fullpix[4];
	get_interp_weights(ring_info_, theta, alpha, fullpix, w);

	pixels = std::vector<long>(4, -1);
	weights = std::vector<double>(4, 0);
	for (size_t i = 0; i < 4; i++) {
		pixels[i] = fullpix[i];
		weights[i] = w[i];
	}
}

G3SkyMapPtr HealpixSkyMap::Rebin(size_t scale, bool norm) const
{
	if (nside_ % scale != 0)
		log_fatal("Map nside must be a multiple of rebinning scale");

	if (scale <= 1)
		return Clone(true);

	HealpixSkyMapPtr out(new HealpixSkyMap(nside_/scale, is_weighted,
	    is_nested_, coord_ref, units, pol_type));

	if (dense_)
		out->ConvertToDense();
	else if (ring_sparse_)
		out->ConvertToRingSparse();
	else if (!indexed_sparse_)
		return out;

	const size_t scale2 = scale * scale;
	const double sqscal = norm ? scale2 : 1.0;

	for (auto i : *this) {
		if (i.second == 0)
			continue;
		long ip = i.first;
		if (!is_nested_)
			ring2nest(nside_, ip, &ip);
		ip /= scale2;
		if (!is_nested_)
			nest2ring(out->nside_, ip, &ip);
		(*out)[ip] += i.second / sqscal;
	}

	return out;
}

void HealpixSkyMap::SetShiftRa(bool shift)
{
	if (shift == shift_ra_)
		return;

	if (!IsRingSparse()) {
		shift_ra_ = shift;
		free(ring_info_);
		ring_info_ = init_map_info(nside_, is_nested_, shift, 1);
		return;
	}

	map_info *ring_info = init_map_info(nside_, is_nested_, shift, 1);
	SparseMapData *ring_sparse = new SparseMapData(ring_info->nring, ring_info->nring);
	long iring, ringpix;
	for (auto i : *this) {
		if (i.second != 0) {
			get_ring_index(ring_info, i.first, &iring, &ringpix);
			(*ring_sparse)(iring, ringpix) = i.second;
		}
	}

	free(ring_info_);
	delete ring_sparse_;
	ring_info_ = ring_info;
	ring_sparse_ = ring_sparse;
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
		PyErr_SetString(PyExc_ValueError,
		    "Cannot set dense to False. Set ringsparse or "
		    "indexedsparse to True to convert from dense.");
		throw boost::python::error_already_set();
	}
}

static void
HealpixSkyMap_setringsparse(HealpixSkyMap &m, bool v)
{
	if (v)
		m.ConvertToRingSparse();
	else {
		PyErr_SetString(PyExc_ValueError,
		    "Cannot set ringsparse to False. Set indexedsparse or "
		    "dense to True to convert from ringsparse.");
		throw boost::python::error_already_set();
	}
}

static void
HealpixSkyMap_setindexedsparse(HealpixSkyMap &m, bool v)
{
	if (v)
		m.ConvertToIndexedSparse();
	else {
		PyErr_SetString(PyExc_ValueError,
		    "Cannot set indexedsparse to False. Set ringsparse or "
		    "dense to True to convert from indexedsparse.");
		throw boost::python::error_already_set();
	}
}

static boost::python::tuple
HealpixSkyMap_nonzeropixels(const HealpixSkyMap &m)
{
	auto i = std::vector<uint64_t>(); // XXX pointers?
	auto d = std::vector<double>();

	m.NonZeroPixels(i, d);

	return boost::python::make_tuple(i, d);
}

static int
HealpixSkyMap_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	namespace bp = boost::python;

	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	bp::handle<> self(bp::borrowed(obj));
	bp::object selfobj(self);
	HealpixSkyMapPtr sm = bp::extract<HealpixSkyMapPtr>(selfobj)();

	sm->ConvertToDense();

	view->obj = obj;
	view->buf = (void*)&(*sm)[0];
	view->len = sm->size() * sizeof(double);
	view->readonly = 0;
	view->itemsize = sizeof(double);
	if (flags & PyBUF_FORMAT)
		view->format = (char *)"d";
	else
		view->format = NULL;

	// XXX: following leaks small amounts of memory!
	view->shape = new Py_ssize_t;
	view->strides = new Py_ssize_t;

	view->ndim = 1;
	view->shape[0] = sm->size();
	view->strides[0] = view->itemsize;

	view->suboffsets = NULL;

	Py_INCREF(obj);

	return 0;
}

static PyBufferProcs healpixskymap_bufferprocs;

#define HEALPIX_SKY_MAP_DOCSTR \
	"HealpixSkyMap is a G3SkyMap with Healpix" \
	"\n\n" \
	"The other meta information is inherited from G3SkyMap that lives in core. \n\n" \
	"If you find that you need numpy functionality from a HealpixSkyMap,\n"\
	" using np.asarray will convert it to a numpy array without copying the data.\n" \
	" any changes to the resulting numpy array will affect the data stored in the\n" \
	" HealpixSkyMap."



G3_SPLIT_SERIALIZABLE_CODE(HealpixSkyMap);

PYBINDINGS("coordinateutils")
{
	using namespace boost::python;

	// Can't use the normal FRAMEOBJECT code since this inherits
	// from an intermediate class. Expanded by hand here.
	object hsm = class_<HealpixSkyMap, bases<G3SkyMap, G3FrameObject>,
	  HealpixSkyMapPtr>(
	  "HealpixSkyMap", HEALPIX_SKY_MAP_DOCSTR, boost::python::no_init)
	    .def(boost::python::init<const HealpixSkyMap &>())
	    .def_pickle(g3frameobject_picklesuite<HealpixSkyMap>())
	    .def(bp::init<size_t, bool, bool, MapCoordReference,
	       G3Timestream::TimestreamUnits, G3SkyMap::MapPolType, bool>(
	         (bp::arg("nside"),
		  bp::args("is_weighted") = true,
		  bp::args("is_nested") = false,
		  bp::arg("coord_ref") = MapCoordReference::Equatorial,
		  bp::arg("units") = G3Timestream::Tcmb,
		  bp::arg("pol_type") = G3SkyMap::None,
		  bp::arg("shift_ra") = false),
	       "Instantiate a HealpixSkyMap with given nside"))
	    .def(bp::init<boost::python::object, bool, bool,
	       MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType>(
		  (bp::arg("data"),
		   bp::args("is_weighted") = true,
		   bp::args("is_nested") = false,
		   bp::arg("coord_ref") = MapCoordReference::Equatorial,
		   bp::arg("units") = G3Timestream::Tcmb,
		   bp::arg("pol_type") = G3SkyMap::None),
	       "Instantiate a Healpix map from existing data. If the data are "
	       "a single numpy array, assumes this is a dense map. Otherwise, "
	       "pass an (indices, data, nside) tuple."))

	    .def(bp::init<const HealpixSkyMap&>(bp::arg("healpix_map")))
	    .def(bp::init<>())
	    .add_property("nside", &HealpixSkyMap::nside)
	    .add_property("nested", &HealpixSkyMap::nested)
	    .add_property("shift_ra", &HealpixSkyMap::IsRaShifted, HealpixSkyMap_setshiftra,
		"True if the ringsparse representation of the map is stored "
		"with the rings centered at ra = 0 deg, rather than ra = 180 deg.")
	    .add_property("dense", &HealpixSkyMap::IsDense, HealpixSkyMap_setdense,
	        "True if the map is stored with all elements, False otherwise. "
	        "If set to True, converts the map to a dense representation." )
	    .add_property("ringsparse", &HealpixSkyMap::IsRingSparse, HealpixSkyMap_setringsparse,
	        "True if the map is stored as a dense 2D region using ring "
	        "ordering (analogous to FlatSkyMap's sparse mode). "
	        "Ring-sparsity is efficient for dense blocks on a ring-ordered "
	        "map (e.g. a continuous sky region), but is inefficient "
	        "otherwise. It applies only to non-nested maps. "
	        "If set to True, converts the map to this representation." )
	    .add_property("indexedsparse", &HealpixSkyMap::IsIndexedSparse, HealpixSkyMap_setindexedsparse,
	        "True if the map is stored as a list of non-zero pixels "
	        "and values. More efficient than ring-sparse for maps with "
	        "holes or very small filling factors. "
	        "If set to True, converts the map to this representation." )

	    .def("nonzero_pixels", &HealpixSkyMap_nonzeropixels,
	        "Returns a list of the indices of the non-zero pixels in the "
	        "map and a list of the values of those non-zero pixels.")
	;
	register_pointer_conversions<HealpixSkyMap>();

	// Add buffer protocol interface
	PyTypeObject *hsmclass = (PyTypeObject *)hsm.ptr();
	healpixskymap_bufferprocs.bf_getbuffer = HealpixSkyMap_getbuffer;
	hsmclass->tp_as_buffer = &healpixskymap_bufferprocs;
#if PY_MAJOR_VERSION < 3
	hsmclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	implicitly_convertible<HealpixSkyMapPtr, G3SkyMapPtr>();
	implicitly_convertible<HealpixSkyMapConstPtr, G3SkyMapConstPtr>();
	implicitly_convertible<HealpixSkyMapPtr, G3SkyMapConstPtr>();
	implicitly_convertible<HealpixSkyMapConstPtr, G3SkyMapConstPtr>();
}

