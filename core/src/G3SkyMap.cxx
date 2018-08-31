#include <pybindings.h>
#include <serialization.h>
#include <G3SkyMap.h>

namespace bp=boost::python;

template <class A> void
G3SkyMap::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("coord_ref", coord_ref);
	ar & cereal::make_nvp("units", units);
	ar & cereal::make_nvp("data", data_);
	ar & cereal::make_nvp("xpix", xpix_);
	ar & cereal::make_nvp("ypix", ypix_);

	ar & cereal::make_nvp("pol_type", pol_type);
	ar & cereal::make_nvp("is_weighted", is_weighted);
}

template <class A> void
G3SkyMapWeights::serialize(A &ar, unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	ar & cereal::make_nvp("TT", TT);
	ar & cereal::make_nvp("TQ", TQ);
	ar & cereal::make_nvp("TU", TU);
	ar & cereal::make_nvp("QQ", QQ);
	ar & cereal::make_nvp("QU", QU);
	ar & cereal::make_nvp("UU", UU);
	if (v > 1) {
		ar & cereal::make_nvp("weight_type", weight_type);
	}
}

G3_SERIALIZABLE_CODE(G3SkyMap);
G3_SERIALIZABLE_CODE(G3SkyMapWeights);

G3SkyMapWeights::G3SkyMapWeights(G3SkyMapConstPtr ref, 
  G3SkyMapWeights::WeightType wt = G3SkyMapWeights::Wpol) :
   TT(ref->Clone(false)), TQ(ref->Clone(false)), TU(ref->Clone(false)),
   QQ(ref->Clone(false)), QU(ref->Clone(false)), UU(ref->Clone(false)),
   weight_type(wt)
{
}

G3SkyMapWeights::G3SkyMapWeights(const G3SkyMapWeights &r) :
    TT(r.TT->Clone(true)), TQ(r.TQ->Clone(true)), TU(r.TU->Clone(true)),
    QQ(r.QQ->Clone(true)), QU(r.QU->Clone(true)), UU(r.UU->Clone(true)),
    weight_type(r.weight_type)
{
}

G3SkyMap::G3SkyMap(bp::object v, MapCoordReference coords, bool is_weighted,
    G3Timestream::TimestreamUnits u, G3SkyMap::MapPolType pol_type) :
    coord_ref(coords), units(u), pol_type(pol_type), is_weighted(is_weighted)
{
	Py_buffer view;
	if (PyObject_GetBuffer(v.ptr(), &view,
	    PyBUF_FORMAT | PyBUF_ANY_CONTIGUOUS) != -1) {
		if (view.ndim == 1){
			xpix_ = view.shape[0];
			ypix_ = 1;
		} else if (view.ndim == 2) {
			ypix_ = view.shape[0];
			xpix_ = view.shape[1];
		} else {
			log_fatal("Only 1 and 2-D maps supported");
		}
		data_.resize(xpix_*ypix_+1);

		if (strcmp(view.format, "d") == 0) {
			memcpy(&(data_[0]), view.buf, view.len);
		} else if (strcmp(view.format, "f") == 0) {
			for (size_t i = 0; i < view.len/sizeof(float); i++)
				data_[i] = ((float *)view.buf)[i];
		} else if (strcmp(view.format, "i") == 0) {
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				data_[i] = ((int *)view.buf)[i];
		} else if (strcmp(view.format, "I") == 0) {
			for (size_t i = 0; i < view.len/sizeof(int); i++)
				data_[i] = ((unsigned int *)view.buf)[i];
		} else if (strcmp(view.format, "l") == 0) {
			for (size_t i = 0; i < view.len/sizeof(long); i++)
				data_[i] = ((unsigned long *)view.buf)[i];
		} else {
			log_fatal("Unknown type code %s", view.format);
		}
		PyBuffer_Release(&view);

		data_[data_.size()-1] = 0;
		return;
	}

	throw bp::error_already_set();
}

size_t G3SkyMap::angle_to_pixel(double alpha, double delta) const {
	std::vector<double> alphas = {alpha};
	std::vector<double> deltas = {delta};
	return angles_to_pixels(alphas, deltas)[0];
}

std::vector<double> G3SkyMap::pixel_to_angle(size_t x_pix, size_t y_pix) const {
	return pixel_to_angle(pixat(x_pix, y_pix));
}

static int
G3SkyMap_getbuffer(PyObject *obj, Py_buffer *view, int flags)
{
	if (view == NULL) {
		PyErr_SetString(PyExc_ValueError, "NULL view");
		return -1;
	}

	view->shape = NULL;

	bp::handle<> self(bp::borrowed(obj));
	bp::object selfobj(self);
	G3SkyMapPtr sm = bp::extract<G3SkyMapPtr>(selfobj)();

	sm->EnsureAllocated();

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

	if (sm->ydim() == 1) {
		view->ndim = 1;
		view->shape[0] = sm->xdim();
		view->strides[0] = view->itemsize;
	} else {
		view->ndim = 2;
		view->shape[0] = sm->ydim();
		view->shape[1] = sm->xdim();
		view->strides[0] = sm->xdim()*view->itemsize;
		view->strides[1] = view->itemsize;
	}

	view->suboffsets = NULL;

	Py_INCREF(obj);

	return 0;
}

static PyBufferProcs skymap_bufferprocs;

static double
skymap_getitem(G3SkyMap &skymap, int i)
{
	skymap.EnsureAllocated();

	if (i < 0)
		i = skymap.size() + i;
	if (size_t(i) >= skymap.size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		bp::throw_error_already_set();
	}

	return skymap[i];
}

static double
skymap_getitem_2d(G3SkyMap &skymap, bp::tuple coords)
{
	int y = bp::extract<int>(coords[0]);
	int x = bp::extract<int>(coords[1]);
	if (x < 0)
		x = skymap.xdim() + x;
	if (y < 0)
		y = skymap.xdim() + y;
	if (size_t(x) >= skymap.xdim()) {
		PyErr_SetString(PyExc_IndexError, "X index out of range");
		bp::throw_error_already_set();
	}
	if (size_t(y) >= skymap.ydim()) {
		PyErr_SetString(PyExc_IndexError, "Y index out of range");
		bp::throw_error_already_set();
	}

	skymap.EnsureAllocated();

	return skymap[skymap.pixat(x, y)];
}

static void
skymap_setitem(G3SkyMap &skymap, int i, double val)
{
	skymap.EnsureAllocated();

	if (i < 0)
		i = skymap.size() + i;
	if (size_t(i) >= skymap.size()) {
		PyErr_SetString(PyExc_IndexError, "Index out of range");
		bp::throw_error_already_set();
	}

	skymap[i] = val;
}

static void
skymap_setitem_2d(G3SkyMap &skymap, bp::tuple coords, double val)
{
	int y = bp::extract<int>(coords[0]);
	int x = bp::extract<int>(coords[1]);
	if (x < 0)
		x = skymap.xdim() + x;
	if (y < 0)
		y = skymap.xdim() + y;

	skymap.EnsureAllocated();

	if (size_t(x) >= skymap.xdim()) {
		PyErr_SetString(PyExc_IndexError, "X index out of range");
		bp::throw_error_already_set();
	}
	if (size_t(y) >= skymap.ydim()) {
		PyErr_SetString(PyExc_IndexError, "Y index out of range");
		bp::throw_error_already_set();
	}

	skymap[skymap.pixat(x, y)] = val;
}

static bp::tuple
skymap_shape(G3SkyMap &skymap)
{
	if (skymap.ydim() == 1)
		return bp::make_tuple(skymap.xdim());
	else
		return bp::make_tuple(skymap.ydim(), skymap.xdim());
}

static G3SkyMapPtr
skymap_copy(G3SkyMap &r)
{
	return r.Clone(true);
}

PYBINDINGS("core") {
	bp::enum_<MapCoordReference>("MapCoordReference")
	    .value("Local", Local)
	    .value("Equatorial", Equatorial)
	    .value("Galactic", Galactic)
	;

	bp::enum_<G3SkyMap::MapPolType>("MapPolType")
	    .value("T", G3SkyMap::T)
	    .value("Q", G3SkyMap::Q)
	    .value("U", G3SkyMap::U)
	    .value("none", G3SkyMap::None) // "None" is reserved in python
	;

	bp::enum_<G3SkyMapWeights::WeightType>("WeightType")
	    .value("Wpol", G3SkyMapWeights::Wpol)
	    .value("Wunpol", G3SkyMapWeights::Wunpol)
	    .value("None", G3SkyMapWeights::None)
	;

	bp::object skymap = EXPORT_FRAMEOBJECT(G3SkyMap, no_init,
	  "Base class for 1- and 2-D skymaps of various projections. Usually "
	  "you want a subclass of this (e.g. FlatSkyMap) rather than using it "
	  "directly.")
	    .def(bp::init<MapCoordReference, size_t, size_t, bool,
	      G3Timestream::TimestreamUnits, G3SkyMap::MapPolType, bool>(
	      (bp::arg("coords"), "xpixels", bp::arg("ypixels")=1,
	       bp::arg("is_weighted")=true,
	       bp::arg("units")=G3Timestream::Tcmb,
	       bp::arg("pol_type")=G3SkyMap::MapPolType::None,
	       bp::arg("allocate_map")=false), "Initialize a sky map with all "
	       "pixels set to zero and the given dimensions"))
	    .def(bp::init<bp::object, MapCoordReference, bool,
	      G3Timestream::TimestreamUnits, G3SkyMap::MapPolType>
	      ((bp::arg("array"), "coords", bp::arg("is_weighted")=true,
	        bp::arg("units")=G3Timestream::Tcmb,
	        bp::arg("pol_type")=G3SkyMap::MapPolType::None),
	      "Create a skymap from a numpy array"))
	    .def_readwrite("coord_ref", &G3SkyMap::coord_ref,
	      "Coordinate system (core.MapCoordReference) of the map (e.g. "
	      "Galactic, Equatorial, etc.)")
	    .def_readwrite("pol_type", &G3SkyMap::pol_type,
	      "Polarization type (core.MapPolType) of the map "
	      "(e.g. core.MapPolType.Q)")
	    .def_readwrite("units", &G3SkyMap::units,
	      "Unit class (core.G3TimestreamUnits) of the map (e.g. "
	      "core.G3TimestreamUnits.Tcmb). Within each unit class, further "
	      "conversions, for example from K to uK, should use core.G3Units.")
	    .def_readwrite("is_weighted", &G3SkyMap::is_weighted,
	      "True if map is multiplied by weights")
	    .add_property("size", &G3SkyMap::size, "Number of pixels in map")
	    .add_property("shape", &skymap_shape, "Shape of map")
	    .add_property("overflow", &G3SkyMap::GetOverflow,
              &G3SkyMap::SetOverflow, "Combined value of data processed by "
	      "the map maker but outside of the map area")
	    .def("__getitem__", &skymap_getitem)
	    .def("__setitem__", &skymap_setitem)
	    .def("__getitem__", &skymap_getitem_2d)
	    .def("__setitem__", &skymap_setitem_2d)
	    .def("__copy__", &skymap_copy)
	    .def("__deepcopy__", &skymap_copy)

	    .def("angles_to_pixels", &G3SkyMap::angles_to_pixels)
	    .def("pixel_to_angle", 
	      (std::vector<double> (G3SkyMap::*)(size_t, size_t) const) 
	      &G3SkyMap::pixel_to_angle, ((bp::arg("x_pix"), bp::arg("y_pix"))))
	    .def("pixel_to_angle", 
	      (std::vector<double> (G3SkyMap::*)(size_t) const) 
	      &G3SkyMap::pixel_to_angle, ((bp::arg("pixel"))))
	    .def("angle_to_pixel", &G3SkyMap::angle_to_pixel,
	      ((bp::arg("alpha"), bp::arg("delta"))))


	;

	// Add buffer protocol interface
	PyTypeObject *smclass = (PyTypeObject *)skymap.ptr();
	skymap_bufferprocs.bf_getbuffer = G3SkyMap_getbuffer;
	smclass->tp_as_buffer = &skymap_bufferprocs;
#if PY_MAJOR_VERSION < 3
	smclass->tp_flags |= Py_TPFLAGS_HAVE_NEWBUFFER;
#endif

	EXPORT_FRAMEOBJECT(G3SkyMapWeights, init<>(), "generic sky weight")
	    .def(bp::init<G3SkyMapConstPtr, G3SkyMapWeights::WeightType>(
	      (bp::arg("skymap"),
               bp::arg("weight_type") = G3SkyMapWeights::Wpol)))
	    .def_readwrite("TT",&G3SkyMapWeights::TT)
	    .def_readwrite("TQ",&G3SkyMapWeights::TQ)
	    .def_readwrite("TU",&G3SkyMapWeights::TU)
	    .def_readwrite("QQ",&G3SkyMapWeights::QQ)
	    .def_readwrite("QU",&G3SkyMapWeights::QU)
	    .def_readwrite("UU",&G3SkyMapWeights::UU)
	    .def_readwrite("weight_type", &G3SkyMapWeights::weight_type)
	    
	    .def("Clone", &G3SkyMapWeights::Clone)
	;
	register_pointer_conversions<G3SkyMapWeights>();


}

