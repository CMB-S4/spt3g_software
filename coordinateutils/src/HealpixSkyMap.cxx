#include <pybindings.h>
#include <serialization.h>
#include <typeinfo>

#include <coordinateutils/HealpixSkyMap.h>
#include <coordinateutils/chealpix.h>

#include "mapdata.h"

HealpixSkyMap::HealpixSkyMap(size_t nside, bool is_weighted,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type) :
      G3SkyMap(coord_ref, is_weighted, u, pol_type),
      nside_(nside), dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
	ring_info_ = init_map_info(nside, 1);
	is_nested_ = false;
}

HealpixSkyMap::HealpixSkyMap(boost::python::object v, bool is_weighted,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type) :
      G3SkyMap(coord_ref, is_weighted, u, pol_type),
      dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
	Py_buffer view;

	is_nested_ = false;

	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
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

		ring_info_ = init_map_info(nside_, 1);

		return;
	}

	throw bp::error_already_set();
}

HealpixSkyMap::HealpixSkyMap(boost::python::object index,
    boost::python::object v, size_t nside, bool is_weighted,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type) :
      G3SkyMap(coord_ref, is_weighted, u, pol_type),
      dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
	Py_buffer indexview, dataview;

	is_nested_ = false;

	if (PyObject_GetBuffer(index.ptr(), &indexview,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1)
		throw bp::error_already_set();

	if (PyObject_GetBuffer(v.ptr(), &dataview,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) == -1) {
		PyBuffer_Release(&indexview);
		throw bp::error_already_set();
	}

	if (indexview.ndim != 1 || dataview.ndim != 1) {
		PyBuffer_Release(&indexview);
		PyBuffer_Release(&dataview);
		log_fatal("Only 1-D maps supported");
	}

	if (indexview.shape[0] != dataview.shape[0]) {
		PyBuffer_Release(&indexview);
		PyBuffer_Release(&dataview);
		log_fatal("Index and data must have matching shapes.");
	}

	if (strcmp(indexview.format, "I") != 0 ||
	    strcmp(indexview.format, "i") != 0) {
		PyBuffer_Release(&indexview);
		PyBuffer_Release(&dataview);
		log_fatal("Indices must be integers.");
	}

	if (strcmp(dataview.format, "d") != 0) {
		PyBuffer_Release(&indexview);
		PyBuffer_Release(&dataview);
		log_fatal("Data must be double-precision (float64).");
	}
	nside_ = npix2nside(dataview.shape[0]);
	ring_info_ = init_map_info(nside_, 1);
	indexed_sparse_ = new std::unordered_map<uint32_t, double>;

	for (size_t i = 0; i < indexview.shape[0]; i++)
		(*indexed_sparse_)[((unsigned int *)indexview.buf)[i]] =
		    ((double *)dataview.buf)[i];
	PyBuffer_Release(&indexview);
	PyBuffer_Release(&dataview);
}

HealpixSkyMap::HealpixSkyMap() :
    G3SkyMap(MapCoordReference::Local, false), nside_(0),
    dense_(NULL), ring_sparse_(NULL), indexed_sparse_(NULL)
{
	ring_info_ = init_map_info(nside_, 1);
	is_nested_ = false;
}

HealpixSkyMap::HealpixSkyMap(const HealpixSkyMap & fm) :
    G3SkyMap(fm), nside_(fm.nside_), dense_(NULL), ring_sparse_(NULL),
    indexed_sparse_(NULL)
{
	if (fm.dense_)
		dense_ = new std::vector<double>(*fm.dense_);
	else if (fm.ring_sparse_)
		ring_sparse_ = new SparseMapData(*fm.ring_sparse_);
	else if (fm.indexed_sparse_)
		indexed_sparse_ = new std::unordered_map<uint32_t, double>(
		    *fm.indexed_sparse_);
	ring_info_ = init_map_info(nside_, 1);
	is_nested_ = false;
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

	free_map_info(ring_info_);
	ring_info_ = init_map_info(nside_, 1);

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
		indexed_sparse_ = new std::unordered_map<uint32_t, double>;
		ar & make_nvp("data", *indexed_sparse_);
		break;
	}
}

G3SkyMapPtr
HealpixSkyMap::Clone(bool copy_data) const
{
	if (copy_data)
		return boost::make_shared<HealpixSkyMap>(*this);
	else
		return boost::make_shared<HealpixSkyMap>(nside_, is_weighted,
		    coord_ref, units, pol_type);
}

double
HealpixSkyMap::operator [] (int i) const
{
	if (dense_)
		return (*dense_)[i];
	if (ring_sparse_) {
		for (long j = 0; j < ring_info_->nring; j++) {
			if (i >= ring_info_->rings[j].startpix &&
			    i < ring_info_->rings[j].startpix +
			    ring_info_->rings[j].ringpix)
				return (*ring_sparse_)(j,
				    ring_info_->rings[j].startpix);
		}
	}
	if (indexed_sparse_)
		return (*indexed_sparse_)[i];
	return 0;
}

double &
HealpixSkyMap::operator [] (int i)
{
	if (dense_)
		return (*dense_)[i];

	if (ring_sparse_) {
		for (long j = 0; j < ring_info_->nring; j++) {
			if (i >= ring_info_->rings[j].startpix &&
			    i < ring_info_->rings[j].startpix +
			    ring_info_->rings[j].ringpix)
				return (*ring_sparse_)(j,
				    ring_info_->rings[j].startpix);
		}
	}

	if (indexed_sparse_)
		return (*indexed_sparse_)[i];

	indexed_sparse_ = new std::unordered_map<uint32_t, double>;
	return (*indexed_sparse_)[i];
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

std::vector<size_t>
HealpixSkyMap::shape() const
{
	std::vector<size_t> shape(1);
	shape[0] = nside2npix(nside_);
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

	if (is_nested_) {
		for (size_t i = 0; i < 4; i++) {
			ring2nest(nside_, fullpix[i], fullpix + i);
		}
	}

	pixels = std::vector<long>(4, -1);
	weights = std::vector<double>(4, 0);
	for (size_t i = 0; i < 4; i++) {
		pixels[i] = fullpix[i];
		weights[i] = w[i];
	}
}

G3SkyMapPtr HealpixSkyMap::rebin(size_t scale) const
{
	// XXX
	return G3SkyMapPtr();
}

#if 0
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

static bool
flatskymap_pysparsity_get(const HealpixSkyMap &fsm)
{
	return !fsm.IsDense();
}

static void
flatskymap_pysparsity_set(HealpixSkyMap &fsm, bool sparse)
{
	if (sparse)
		fsm.ConvertToSparse();
	else
		fsm.ConvertToDense();
}
#endif

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
	object fsm = class_<HealpixSkyMap, bases<G3SkyMap>, HealpixSkyMapPtr>(
	  "HealpixSkyMap", HEALPIX_SKY_MAP_DOCSTR, boost::python::no_init)
	    .def(boost::python::init<const HealpixSkyMap &>())
	    .def_pickle(g3frameobject_picklesuite<HealpixSkyMap>())
	    .def(bp::init<size_t, bool, MapCoordReference,
	       G3Timestream::TimestreamUnits, G3SkyMap::MapPolType>(
	         (bp::arg("nside"),
		  bp::args("is_weighted") = true,
		  bp::arg("coord_ref") = MapCoordReference::Equatorial,
		  bp::arg("units") = G3Timestream::Tcmb,
		  bp::arg("pol_type") = G3SkyMap::None)))
#if 0
	    .def(bp::init<boost::python::object, bool,
	       MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType>(
		  (bp::arg("dense_healpix"),
	           bp::args("is_weighted") = true,
		   bp::arg("coord_ref") = MapCoordReference::Equatorial,
		   bp::arg("units") = G3Timestream::Tcmb,
		   bp::arg("pol_type") = G3SkyMap::None)))
	    .def(bp::init<boost::python::object, boost::python::object, size_t,
	       bool, MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType>(
		  (bp::arg("indices"), bp::arg("data"), bp::arg("nside"),
	           bp::args("is_weighted") = true,
		   bp::arg("coord_ref") = MapCoordReference::Equatorial,
		   bp::arg("units") = G3Timestream::Tcmb,
		   bp::arg("pol_type") = G3SkyMap::None)))
#endif

	    .def(bp::init<const HealpixSkyMap&>(bp::arg("healpix_map")))
	    .def(bp::init<>())
#if 0
	    .add_property("sparse", flatskymap_pysparsity_get, flatskymap_pysparsity_set,
	       "True if the map is stored with column and row zero-suppression, False if "
	       "every pixel is stored. Map sparsity can be changed by setting this to True "
	       "(or False).")
#endif

	;
	register_pointer_conversions<HealpixSkyMap>();

#if 0
	// Add buffer protocol interface
	PyTypeObject *fsmclass = (PyTypeObject *)fsm.ptr();
	flatskymap_bufferprocs.bf_getbuffer = HealpixSkyMap_getbuffer;
	fsmclass->tp_as_buffer = &flatskymap_bufferprocs;
#if PY_MAJOR_VERSION < 3
	fsmclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif
#endif

	implicitly_convertible<HealpixSkyMapPtr, G3SkyMapPtr>();
	implicitly_convertible<HealpixSkyMapConstPtr, G3SkyMapConstPtr>();
	implicitly_convertible<HealpixSkyMapPtr, G3SkyMapConstPtr>();
	implicitly_convertible<HealpixSkyMapConstPtr, G3SkyMapConstPtr>();
}

