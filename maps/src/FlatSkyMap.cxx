#include <pybindings.h>
#include <serialization.h>
#include <typeinfo>

#include <maps/FlatSkyMap.h>
#include <maps/FlatSkyProjection.h>

#include "mapdata.h"

FlatSkyMap::FlatSkyMap(size_t x_len, size_t y_len, double res, bool weighted,
    MapProjection proj, double alpha_center, double delta_center,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, double x_res, double x_center, double y_center,
    bool flat_pol, G3SkyMap::MapPolConv pol_conv) :
      G3SkyMap(coord_ref, weighted, u, pol_type, pol_conv),
      proj_info(x_len, y_len, res, alpha_center, delta_center, x_res, proj, x_center, y_center),
      dense_(NULL), sparse_(NULL), xpix_(x_len), ypix_(y_len), flat_pol_(flat_pol)
{
}

FlatSkyMap::FlatSkyMap(boost::python::object v, double res,
    bool weighted, MapProjection proj,
    double alpha_center, double delta_center,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, double x_res, double x_center, double y_center,
    bool flat_pol, G3SkyMap::MapPolConv pol_conv) :
      G3SkyMap(coord_ref, weighted, u, pol_type, pol_conv),
      dense_(NULL), sparse_(NULL), flat_pol_(flat_pol)
{
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_C_CONTIGUOUS) != -1) {
		if (view.ndim == 2) {
			ypix_ = view.shape[0];
			xpix_ = view.shape[1];
		} else {
			log_fatal("Only 2-D maps supported");
		}
		proj_info = FlatSkyProjection(xpix_, ypix_, res, alpha_center,
		    delta_center, x_res, proj, x_center, y_center);
		ConvertToDense();

		double *d = &(*dense_)(0,0);
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

		return;
	}

	throw bp::error_already_set();
}

FlatSkyMap::FlatSkyMap(const FlatSkyProjection & fp,
    MapCoordReference coord_ref, bool weighted,
    G3Timestream::TimestreamUnits u, G3SkyMap::MapPolType pol_type,
    bool flat_pol, G3SkyMap::MapPolConv pol_conv) :
      G3SkyMap(coord_ref, weighted, u, pol_type, pol_conv),
      proj_info(fp), dense_(NULL), sparse_(NULL),
      xpix_(fp.xdim()), ypix_(fp.ydim()), flat_pol_(flat_pol)
{
}

FlatSkyMap::FlatSkyMap() :
    G3SkyMap(MapCoordReference::Local, 0), proj_info(0, 0, 0),
    dense_(NULL), sparse_(NULL), xpix_(0), ypix_(0), flat_pol_(false)
{
}

FlatSkyMap::FlatSkyMap(const FlatSkyMap & fm) :
    G3SkyMap(fm), proj_info(fm.proj_info), dense_(NULL), sparse_(NULL),
    xpix_(fm.xpix_), ypix_(fm.ypix_), flat_pol_(fm.flat_pol_)
{
	if (fm.dense_)
		dense_ = new DenseMapData(*fm.dense_);
	else if (fm.sparse_)
		sparse_ = new SparseMapData(*fm.sparse_);
}

FlatSkyMap::~FlatSkyMap()
{
	if (dense_)
		delete dense_;
	if (sparse_)
		delete sparse_;
}

template <class A> void
FlatSkyMap::save(A &ar, unsigned v) const
{
	using namespace cereal;

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
	ar & make_nvp("proj_info", proj_info);
	ar & make_nvp("xpix", xpix_);
	ar & make_nvp("ypix", ypix_);
	if (dense_) {
		ar & make_nvp("store", 2);
		ar & make_nvp("data", *dense_);
	} else if (sparse_) {
		ar & make_nvp("store", 1);
		ar & make_nvp("data", *sparse_);
	} else {
		ar & make_nvp("store", 0);
	}
	ar & make_nvp("flat_pol", flat_pol_);
}

template <class A> void
FlatSkyMap::load(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));

	if (v >= 2) {
		ar & make_nvp("proj_info", proj_info);
	} else {
		MapProjection proj;
		double alpha_center, delta_center, res, x_res;
		ar & make_nvp("proj", proj);
		ar & make_nvp("alpha_center", alpha_center);
		ar & make_nvp("delta_center", delta_center);
		ar & make_nvp("res", res);
		ar & make_nvp("x_res", x_res);
		proj_info.initialize(xpix_, ypix_, res, alpha_center,
		    delta_center, x_res, proj);
	}

	if (v >= 3) {
		int store;
		ar & make_nvp("xpix", xpix_);
		ar & make_nvp("ypix", ypix_);
		ar & make_nvp("store", store);
		if (dense_) {
			delete dense_;
			dense_ = NULL;
		}
		if (sparse_) {
			delete sparse_;
			sparse_ = NULL;
		}
		switch (store) {
		case 2:
			dense_ = new DenseMapData(0, 0);
			ar & make_nvp("dense", *dense_);
			assert(dense_->xdim() == xpix_);
			assert(dense_->ydim() == ypix_);
			break;
		case 1:
			sparse_ = new SparseMapData(0, 0);
			ar & make_nvp("sparse", *sparse_);
			assert(sparse_->xdim() == xpix_);
			assert(sparse_->ydim() == ypix_);
			break;
		}
	}

	if (v >= 4) {
		ar & make_nvp("flat_pol", flat_pol_);
	} else {
		flat_pol_ = false;
	}
}

void
FlatSkyMap::InitFromV1Data(std::vector<size_t> dims, const std::vector<double> &data)
{
	xpix_ = dims[0];
	ypix_ = dims[1];

	if (data.size() > 0) {
		dense_ = new DenseMapData(dims[0], dims[1]);
		(*dense_) = data;
	}
}

FlatSkyMap::const_iterator::const_iterator(const FlatSkyMap &map, bool begin) :
    map_(map)
{
	if (map_.dense_) {
		auto iter = begin ? map_.dense_->begin() : map_.dense_->end();
		x_ = iter.x;
		y_ = iter.y;
	} else if (map_.sparse_) {
		auto iter = begin ? map_.sparse_->begin() : map_.sparse_->end();
		x_ = iter.x;
		y_ = iter.y;
	} else {
		x_ = 0;
		y_ = 0;
	}

	set_value();
}

FlatSkyMap::const_iterator
FlatSkyMap::const_iterator::operator++()
{
	if (map_.dense_) {
		DenseMapData::const_iterator iter(*map_.dense_, x_, y_);
		++iter;
		x_ = iter.x;
		y_ = iter.y;
	} else if (map_.sparse_) {
		SparseMapData::const_iterator iter(*map_.sparse_, x_, y_);
		++iter;
		x_ = iter.x;
		y_ = iter.y;
	}

	set_value();
	return *this;
}

void
FlatSkyMap::ConvertToDense()
{
	if (dense_)
		return;

	if (sparse_) {
		dense_ = sparse_->to_dense();
		delete sparse_;
		sparse_ = NULL;
	} else {
		// Unallocated
		dense_ = new DenseMapData(xpix_, ypix_);
	}
}

void
FlatSkyMap::ConvertToSparse()
{
	// Do nothing if unallocated (maximumally sparse) or already sparse.
	if (!dense_)
		return;

	sparse_ = new SparseMapData(*dense_);
	delete dense_;
	dense_ = NULL;
}

void FlatSkyMap::Compress(bool zero_nans)
{
	if (dense_) {
		if (zero_nans) {
			for (auto i =  dense_->begin(); i != dense_->end(); i++) {
				if ((*i) != (*i))
					(*dense_)(i.x, i.y) = 0;
			}
		}

		ConvertToSparse();
		return;
	}

	if (!sparse_)
		return;

	auto *old_sparse = sparse_;
	sparse_ = new SparseMapData(xpix_, ypix_);
	for (auto i = old_sparse->begin(); i != old_sparse->end(); i++) {
		if ((*i) == 0 || (zero_nans && (*i) != (*i)))
			continue;
		(*sparse_)(i.x, i.y) = (*i);
	}

	delete old_sparse;
}

G3SkyMapPtr
FlatSkyMap::Clone(bool copy_data) const
{
	if (copy_data)
		return boost::make_shared<FlatSkyMap>(*this);
	else
		return boost::make_shared<FlatSkyMap>(proj_info,
		    coord_ref, weighted, units, pol_type, flat_pol_, pol_conv);
}

double
FlatSkyMap::at(size_t x, size_t y) const
{
	if (dense_)
		return dense_->at(x, y);
	if (sparse_)
		return sparse_->at(x, y);
	return 0;
}

double &
FlatSkyMap::operator () (size_t x, size_t y)
{
	g3_assert(!(x < 0 || x >= xpix_ || y < 0 || y >= ypix_));

	if (dense_)
		return (*dense_)(x, y);
	if (!sparse_)
		sparse_ = new SparseMapData(xpix_, ypix_);
	return (*sparse_)(x, y);
}

double
FlatSkyMap::at(size_t i) const
{
	return this->at(i % xpix_, i / xpix_);
}

double &
FlatSkyMap::operator [] (size_t i)
{
	return (*this)(i % xpix_, i / xpix_);
}

#define flatskymap_arithmetic(op, sparsenull) \
G3SkyMap &FlatSkyMap::operator op(const G3SkyMap &rhs) {\
	g3_assert(IsCompatible(rhs)); \
	try { \
		const FlatSkyMap& b = dynamic_cast<const FlatSkyMap &>(rhs); \
		if (dense_) { \
			if (b.dense_) \
				(*dense_) op (*b.dense_); \
			else if (b.sparse_) \
				(*dense_) op (*b.sparse_); \
			else \
				sparsenull \
		} else if (sparse_) { \
			if (b.dense_) \
				(*sparse_) op (*b.dense_); \
			else if (b.sparse_) \
				(*sparse_) op (*b.sparse_); \
			else \
				sparsenull \
		} else { \
			if (b.dense_) { \
				ConvertToDense(); \
				(*dense_) op (*b.dense_); \
			} else if (b.sparse_) { \
				sparse_ = new SparseMapData(xpix_, ypix_); \
				(*sparse_) op (*b.sparse_); \
			} else { \
				sparsenull \
			} \
		} \
		return *this; \
	} catch (const std::bad_cast& e) { \
		return G3SkyMap::operator op(rhs); \
	} \
}

flatskymap_arithmetic(+=, {})
flatskymap_arithmetic(-=, {})
flatskymap_arithmetic(/=, {ConvertToDense(); (*this->dense_) /= 0.0;})

G3SkyMap &FlatSkyMap::operator *=(const G3SkyMap &rhs) {
	g3_assert(IsCompatible(rhs));
	try {
		const FlatSkyMap& b = dynamic_cast<const FlatSkyMap &>(rhs);
		bool zero = false;
		if (dense_) {
			if (b.dense_)
				(*dense_) *= (*b.dense_);
			else if (b.sparse_)
				(*dense_) *= (*b.sparse_);
			else 
				zero = true;
		} else if (sparse_) {
			if (b.dense_)
				(*sparse_) *= (*b.dense_);
			else if (b.sparse_)
				(*sparse_) *= (*b.sparse_);
			else
				zero = true;
		} else {
			zero = true;
		}

		if (zero) {
			if (sparse_)
				delete sparse_;
			if (dense_)
				delete dense_;
			dense_=NULL;
			sparse_=NULL;
		}

		return *this;
	} catch (const std::bad_cast& e) {
		return G3SkyMap::operator *=(rhs);
	}
}

G3SkyMap &
FlatSkyMap::operator +=(double b)
{
	if (b == 0)
		return *this;

	if (!dense_)
		ConvertToDense();
	(*dense_) += b;
	return *this;
}

G3SkyMap &
FlatSkyMap::operator -=(double b)
{
	if (b == 0)
		return *this;

	if (!dense_)
		ConvertToDense();
	(*dense_) -= b;
	return *this;
}

G3SkyMap &
FlatSkyMap::operator *=(double b)
{
	if (b == 0) {
		if (sparse_)
			delete sparse_;
		if (dense_)
			delete dense_;
		sparse_ = NULL;
		dense_ = NULL;
		return *this;
	}

	if (dense_)
		(*dense_) *= b;
	else if (sparse_)
		(*sparse_) *= b;
	return *this;
}

G3SkyMap &
FlatSkyMap::operator /=(double b)
{
	if (b == 0)
		ConvertToDense();

	if (dense_)
		(*dense_) /= b;
	else if (sparse_)
		(*sparse_) /= b;

	return *this;
}

std::string
FlatSkyMap::Description() const
{
	std::ostringstream os;

	os.precision(1);

	os << proj_info.Description() << " in ";

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

bool FlatSkyMap::IsCompatible(const G3SkyMap & other) const {
	try {
		const FlatSkyMap & flat = dynamic_cast<const FlatSkyMap&>(other);
		return (G3SkyMap::IsCompatible(other) &&
			proj_info.IsCompatible(flat.proj_info));
	} catch(const std::bad_cast& e) {
		return false;
	}
}

void
FlatSkyMap::NonZeroPixels(std::vector<uint64_t> &indices,
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

std::vector<size_t> FlatSkyMap::shape() const {
	return {xpix_, ypix_};
}

size_t FlatSkyMap::NpixAllocated() const {
	if (dense_)
		return xpix_*ypix_;
	if (sparse_)
		return sparse_->nonzero();
	return 0;
}

#define GETSET(name, cname, type)              \
	type FlatSkyMap::name() const          \
	{                                      \
		return proj_info.name();       \
	}                                      \
	void FlatSkyMap::Set##cname(type name) \
	{                                      \
		proj_info.Set##cname(name);    \
	}

GETSET(proj, Proj, MapProjection);
GETSET(alpha_center, AlphaCenter, double);
GETSET(delta_center, DeltaCenter, double);
GETSET(x_center, XCenter, double);
GETSET(y_center, YCenter, double);
GETSET(xres, XRes, double);
GETSET(yres, YRes, double);
GETSET(res, Res, double);

std::vector<double> FlatSkyMap::AngleToXY(double alpha, double delta) const {
	return proj_info.AngleToXY(alpha, delta);
}

std::vector<double> FlatSkyMap::XYToAngle(double x, double y) const {
	return proj_info.XYToAngle(x, y, false);
}

size_t FlatSkyMap::AngleToPixel(double alpha, double delta) const {
	return proj_info.AngleToPixel(alpha, delta);
}

std::vector<double> FlatSkyMap::PixelToAngle(size_t pixel) const {
	return proj_info.PixelToAngle(pixel, false);
}

std::vector<double> FlatSkyMap::PixelToAngle(size_t x, size_t y) const {
	return proj_info.PixelToAngle(y*xpix_ + x, false);
}

std::vector<double> FlatSkyMap::PixelToAngleWrapRa(size_t pixel) const {
	return proj_info.PixelToAngle(pixel, true);
}

std::vector<double> FlatSkyMap::PixelToAngleGrad(size_t pixel, double h) const {
	return proj_info.PixelToAngleGrad(pixel, h);
}

void FlatSkyMap::GetRebinAngles(long pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas) const
{
	proj_info.GetRebinAngles(pixel, scale, alphas, deltas, false);
}

void FlatSkyMap::GetInterpPixelsWeights(double alpha, double delta,
    std::vector<long> & pixels, std::vector<double> & weights) const
{
	proj_info.GetInterpPixelsWeights(alpha, delta, pixels, weights);
}

G3SkyMapPtr FlatSkyMap::Rebin(size_t scale, bool norm) const
{
	const double sqscal = norm ? scale*scale : 1.0;
	if ((xpix_ % scale != 0) || (ypix_ % scale != 0)) {
		log_fatal("Map dimensions must be a multiple of rebinning scale");
	}

	if (scale <= 1)
		return Clone(true);

	FlatSkyProjection p(proj_info.Rebin(scale));
	FlatSkyMapPtr out(new FlatSkyMap(p, coord_ref, weighted, units, pol_type,
	    flat_pol_, pol_conv));

	if (dense_)
		out->ConvertToDense();
	else if (!sparse_)
		return out;

	for (auto i : *this) {
		if (i.second == 0)
			continue;
		size_t x = (i.first % xpix_) / scale;
		size_t y = ((size_t)(i.first / xpix_)) / scale;
		size_t ip = x + y * out->xpix_;
		(*out)[ip] += i.second / sqscal;
	}

	return out;
}

G3SkyMapPtr FlatSkyMap::ExtractPatch(size_t x0, size_t y0, size_t width, size_t height) const
{
	g3_assert(x0 >= 0 && x0 + width <= xpix_);
	g3_assert(y0 >= 0 && y0 + height <= ypix_);

	FlatSkyProjection p(proj_info.OverlayPatch(x0, y0, width, height));
	FlatSkyMapPtr out(new FlatSkyMap(p, coord_ref, weighted, units, pol_type,
	    flat_pol_, pol_conv));

	for (size_t x = 0; x < width; x++) {
		for (size_t y = 0; y < height; y++) {
			double v = this->at(x + x0, y + y0);
			if (v != 0)
				(*out)(x, y) = v;
		}
	}

	return out;
}

void FlatSkyMap::InsertPatch(const FlatSkyMap &patch)
{
	std::vector<double> loc = proj_info.GetPatchLocation(patch.proj_info);
	double x0 = loc[0];
	double y0 = loc[1];
	g3_assert(x0 >= 0 && x0 + patch.xpix_ <= xpix_);
	g3_assert(y0 >= 0 && y0 + patch.ypix_ <= ypix_);

	for (auto i : patch) {
		size_t x = (size_t)(i.first % patch.xpix_) + x0;
		size_t y = (size_t)(i.first / patch.xpix_) + y0;
		(*this)(x, y) = i.second;
	}
}

G3SkyMapPtr FlatSkyMap::Pad(size_t width, size_t height, double fill) const
{
	g3_assert(width > xpix_);
	g3_assert(height > ypix_);

	double x0 = (width - xpix_) / 2;
	double y0 = (height - ypix_) / 2;
	FlatSkyProjection p(proj_info.OverlayPatch(x0, y0, width, height));
	FlatSkyMapPtr out(new FlatSkyMap(p, coord_ref, weighted, units, pol_type,
	    flat_pol_, pol_conv));

	if (fill != 0)
		(*out) += fill;

	for (auto i : *this) {
		size_t x = (size_t)(i.first % xpix_) + x0;
		size_t y = (size_t)(i.first / xpix_) + y0;
		if (i.second != fill)
			(*out)(x, y) = i.second;
	}

	return out;
}

static int
FlatSkyMap_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	namespace bp = boost::python;

	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	bp::handle<> self(bp::borrowed(obj));
	bp::object selfobj(self);
	FlatSkyMapPtr sm = bp::extract<FlatSkyMapPtr>(selfobj)();

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
	view->shape = new Py_ssize_t[2];
	view->strides = new Py_ssize_t[2];

	view->ndim = 2;
	view->shape[0] = sm->shape()[1]; // Numpy has swapped indexing
	view->shape[1] = sm->shape()[0];
	view->strides[0] = sm->shape()[0]*view->itemsize;
	view->strides[1] = view->itemsize;

	view->suboffsets = NULL;

	Py_INCREF(obj);

	return 0;
}

static PyBufferProcs flatskymap_bufferprocs;

static double
flatskymap_getitem_2d(const FlatSkyMap &skymap, bp::tuple coords)
{
	ssize_t y = bp::extract<ssize_t>(coords[0]); // Swapped to match numpy
	ssize_t x = bp::extract<ssize_t>(coords[1]);
	if (x < 0)
		x = skymap.shape()[0] + x;
	if (y < 0)
		y = skymap.shape()[1] + y;
	if (size_t(x) >= skymap.shape()[0]) {
		PyErr_SetString(PyExc_IndexError, "X index out of range");
		bp::throw_error_already_set();
	}
	if (size_t(y) >= skymap.shape()[1]) {
		PyErr_SetString(PyExc_IndexError, "Y index out of range");
		bp::throw_error_already_set();
	}

	return skymap.at(x, y);
}

static double
flatskymap_setitem_2d(FlatSkyMap &skymap, bp::tuple coords, double val)
{
	ssize_t y = bp::extract<ssize_t>(coords[0]); // Swapped to match numpy
	ssize_t x = bp::extract<ssize_t>(coords[1]);
	if (x < 0)
		x = skymap.shape()[0] + x;
	if (y < 0)
		y = skymap.shape()[1] + y;
	if (size_t(x) >= skymap.shape()[0]) {
		PyErr_SetString(PyExc_IndexError, "X index out of range");
		bp::throw_error_already_set();
	}
	if (size_t(y) >= skymap.shape()[1]) {
		PyErr_SetString(PyExc_IndexError, "Y index out of range");
		bp::throw_error_already_set();
	}

	skymap(x, y) = val;
	return val;
}

static double
flatskymap_getitem_1d(const G3SkyMap &skymap, size_t i)
{

	if (i < 0)
		i = skymap.size() + i;
	if (size_t(i) >= skymap.size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		bp::throw_error_already_set();
	}

	return skymap.at(i);
}

static void
flatskymap_setitem_1d(G3SkyMap &skymap, size_t i, double val)
{

	if (i < 0)
		i = skymap.size() + i;
	if (size_t(i) >= skymap.size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		bp::throw_error_already_set();
	}

	skymap[i] = val;
}

static bool
flatskymap_pysparsity_get(const FlatSkyMap &fsm)
{
	return !fsm.IsDense();
}

static void
flatskymap_pysparsity_set(FlatSkyMap &fsm, bool sparse)
{
	if (sparse)
		fsm.ConvertToSparse();
	else
		fsm.ConvertToDense();
}

static boost::python::tuple
flatskymap_nonzeropixels(const FlatSkyMap &m)
{
	auto i = std::vector<uint64_t>(); // XXX pointers?
	auto d = std::vector<double>();

	m.NonZeroPixels(i, d);

	return boost::python::make_tuple(i, d);
}


G3_SPLIT_SERIALIZABLE_CODE(FlatSkyMap);

PYBINDINGS("maps")
{
	using namespace boost::python;

	// Can't use the normal FRAMEOBJECT code since this inherits
	// from an intermediate class. Expanded by hand here.
	object fsm = class_<FlatSkyMap, bases<G3SkyMap, G3FrameObject>,
	  FlatSkyMapPtr>("FlatSkyMap",
	  "FlatSkyMap is a G3SkyMap with the extra meta information about the "
	  "particular flat sky projection included.  In practice it behaves "
	  "(mostly) like a 2d numpy array.  The pixels are normally indexed with "
	  "an 1d pixel index.  If you find that you need numpy functionality from "
	  "a FlatSkyMap, e.g. for slicing across the two dimensions, you can "
	  "access a numpy representation of the map using `np.asarray(m)`. This "
	  "does not copy the data, so any changes to the resulting array will "
	  "affect the data stored in the map.", boost::python::no_init)
	    .def(boost::python::init<const FlatSkyMap &>())
	    .def_pickle(g3frameobject_picklesuite<FlatSkyMap>())
	    .def(bp::init<size_t, size_t, double, bool, MapProjection, double,
	       double, MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType, double, double, double, bool,
	       G3SkyMap::MapPolConv>(
		 (bp::arg("x_len"), bp::arg("y_len"), bp::arg("res"),
		  bp::arg("weighted") = true,
		  bp::arg("proj") = MapProjection::ProjNone,
		  bp::arg("alpha_center") = 0, bp::arg("delta_center") = 0,
		  bp::arg("coord_ref") = MapCoordReference::Equatorial,
		  bp::arg("units") = G3Timestream::Tcmb,
		  bp::arg("pol_type") = G3SkyMap::None, bp::arg("x_res") = 0,
		  bp::arg("x_center") = 0.0 / 0.0, bp::arg("y_center") = 0.0 / 0.0,
		  bp::arg("flat_pol") = false,
		  bp::arg("pol_conv") = G3SkyMap::ConvNone)))
	    .def(bp::init<boost::python::object, double, bool, MapProjection,
	       double, double, MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType, double, double, double, bool,
	       G3SkyMap::MapPolConv>(
		  (bp::arg("obj"), bp::arg("res"),
		   bp::arg("weighted") = true,
		   bp::arg("proj") = MapProjection::ProjNone,
		   bp::arg("alpha_center") = 0, bp::arg("delta_center") = 0,
		   bp::arg("coord_ref") = MapCoordReference::Equatorial,
		   bp::arg("units") = G3Timestream::Tcmb,
		   bp::arg("pol_type") = G3SkyMap::None,
		   bp::arg("x_res") = 0,
		   bp::arg("x_center") = 0.0 / 0.0, bp::arg("y_center") = 0.0 / 0.0,
		   bp::arg("flat_pol") = false,
		   bp::arg("pol_conv") = G3SkyMap::ConvNone)))

	    .def(bp::init<const FlatSkyMap&>(bp::arg("flat_map")))
	    .def(bp::init<>())
	    .add_property("proj", &FlatSkyMap::proj, &FlatSkyMap::SetProj,
	      "Map projection (one of maps.MapProjection)")
	    .add_property("alpha_center", &FlatSkyMap::alpha_center,
	      &FlatSkyMap::SetAlphaCenter, "Horizontal axis center position")
	    .add_property("delta_center", &FlatSkyMap::delta_center,
	      &FlatSkyMap::SetDeltaCenter, "Vertical axis center position")
	    .add_property("x_center", &FlatSkyMap::x_center,
	      &FlatSkyMap::SetXCenter, "Horizontal axis center pixel position")
	    .add_property("y_center", &FlatSkyMap::y_center,
	      &FlatSkyMap::SetYCenter, "Vertical axis center pixel position")
	    .add_property("res", &FlatSkyMap::res, &FlatSkyMap::SetRes,
	      "Map resolution in angular units for maps with square pixels")
	    .add_property("x_res", &FlatSkyMap::xres, &FlatSkyMap::SetXRes,
	      "Resolution in X direction for maps with rectangular pixels")
	    .add_property("y_res", &FlatSkyMap::yres, &FlatSkyMap::SetYRes,
	      "Resolution in Y direction for maps with rectangular pixels")

	    .def("pixel_to_angle",
	      (std::vector<double> (FlatSkyMap::*)(size_t) const)
	      &FlatSkyMap::PixelToAngle, bp::arg("pixel"),
	      "Compute the sky coordinates of the given 1D pixel")
	    .def("pixel_to_angle",
	      (std::vector<double> (FlatSkyMap::*)(size_t, size_t) const)
	      &FlatSkyMap::PixelToAngle, (bp::arg("x"), bp::arg("y")),
	      "Compute the sky coordinates of the given 2D pixel (also see "
	      "xy_to_angle()")

	    .def("xy_to_angle", &FlatSkyMap::XYToAngle,
	      (bp::arg("x"), bp::arg("y")),
	       "Compute the sky coordinates of the input flat 2D coordinates")
	    .def("angle_to_xy", &FlatSkyMap::AngleToXY,
	      (bp::arg("alpha"), bp::arg("delta")),
	       "Compute the flat 2D coordinates of the input sky coordinates")

	    .add_property("sparse", flatskymap_pysparsity_get, flatskymap_pysparsity_set,
	       "True if the map is stored with column and row zero-suppression, False if "
	       "every pixel is stored. Map sparsity can be changed by setting this to True "
	       "(or False).")

	    .def("nonzero_pixels", &flatskymap_nonzeropixels,
		"Returns a list of the indices of the non-zero pixels in the "
		"map and a list of the values of those non-zero pixels.")

	    .def("extract_patch", &FlatSkyMap::ExtractPatch,
		(bp::arg("x0"), bp::arg("y0"), bp::arg("width"), bp::arg("height")),
		"Returns a map of shape (width, height) containing a rectangular patch "
		"of the parent map.  The (0,0) pixel in the output map corresponds to "
		"(x0, y0) in the parent map, and the angular location of each pixel on "
		"the sky is maintained.")

	    .def("insert_patch", &FlatSkyMap::InsertPatch, (bp::arg("patch")),
		"Inserts a patch (e.g. as extracted using extract_patch) into the "
		"parent map.  The coordinate system and angular center of the patch "
		"must match that of the parent map.")

	    .def("pad", &FlatSkyMap::Pad,
		(bp::arg("width"), bp::arg("height"), bp::arg("fill") = 0),
		"Returns a map of shape (width, height) containing the parent map "
		"centered within it.  The angular location of each pixel on the sky "
		"is maintained.  Empty pixels can be optionally filled.")

	    .add_property("flat_pol", &FlatSkyMap::IsPolFlat, &FlatSkyMap::SetFlatPol,
		"True if this map has been flattened using flatten_pol.")

	    .def("__getitem__", flatskymap_getitem_1d)
	    .def("__setitem__", flatskymap_setitem_1d)
	    .def("__getitem__", flatskymap_getitem_2d)
	    .def("__setitem__", flatskymap_setitem_2d)
	;
	register_pointer_conversions<FlatSkyMap>();

	// Add buffer protocol interface
	PyTypeObject *fsmclass = (PyTypeObject *)fsm.ptr();
	flatskymap_bufferprocs.bf_getbuffer = FlatSkyMap_getbuffer;
	fsmclass->tp_as_buffer = &flatskymap_bufferprocs;
#if PY_MAJOR_VERSION < 3
	fsmclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	implicitly_convertible<FlatSkyMapPtr, G3SkyMapPtr>();
	implicitly_convertible<FlatSkyMapConstPtr, G3SkyMapConstPtr>();
	implicitly_convertible<FlatSkyMapPtr, G3SkyMapConstPtr>();
	implicitly_convertible<FlatSkyMapConstPtr, G3SkyMapConstPtr>();
}

