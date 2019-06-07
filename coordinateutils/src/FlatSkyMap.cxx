#include <pybindings.h>
#include <serialization.h>
#include <typeinfo>

#include <coordinateutils/FlatSkyMap.h>
#include <coordinateutils/flatskyprojection.h>

#include "mapdata.h"

FlatSkyMap::FlatSkyMap(int x_len, int y_len, double res, bool is_weighted,
    MapProjection proj, double alpha_center, double delta_center,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, double x_res) :
      G3SkyMap(coord_ref, x_len, y_len, is_weighted, u, pol_type),
      proj_info(x_len, y_len, res, alpha_center, delta_center, x_res, proj),
      dense_(NULL), sparse_(NULL)
{
}

#if 0
// Needs to pass-through to G3SkyMap constructor, so can't be an
// out-of-class make_constructor() thing
FlatSkyMap::FlatSkyMap(boost::python::object v, double res,
    bool is_weighted, MapProjection proj,
    double alpha_center, double delta_center,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, double x_res) :
      G3SkyMap(v, coord_ref, is_weighted, u, pol_type),
      proj_info(xpix_, ypix_, res, alpha_center, delta_center, x_res, proj), 
      dense_(NULL), sparse_(NULL)
{
}
#endif

FlatSkyMap::FlatSkyMap(const FlatSkyProjection & fp,
    MapCoordReference coord_ref, bool is_weighted,
    G3Timestream::TimestreamUnits u, G3SkyMap::MapPolType pol_type) :
      G3SkyMap(coord_ref, fp.xdim(), fp.ydim(), is_weighted, u, pol_type),
      proj_info(fp)
{
}

FlatSkyMap::FlatSkyMap() :
    G3SkyMap(MapCoordReference::Local, 0), proj_info(0, 0, 0)
{
}

FlatSkyMap::FlatSkyMap(const FlatSkyMap & fm) :
    G3SkyMap(fm), proj_info(fm.proj_info)
{
}

template <class A> void
FlatSkyMap::save(A &ar, unsigned v) const
{
	using namespace cereal;

	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
	ar & make_nvp("proj_info", proj_info);
	if (dense_) {
		ar & make_nvp("store", 2);
		ar & make_nvp("data", *dense_);
	} else if (sparse_) {
		ar & make_nvp("store", 1);
		ar & make_nvp("data", *sparse_);
	} else {
		ar & make_nvp("store", 0);
	}
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
		proj_info.initialize(xpix_, ypix_, res, alpha_center, delta_center, x_res, proj);
	}

	if (v >= 3) {
		int store;
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
			dense_ = new DenseMapData(xpix_, ypix_);
			ar & make_nvp("dense", *dense_);
			break;
		case 1:
			sparse_ = new SparseMapData(xpix_, ypix_);
			ar & make_nvp("sparse", *sparse_);
			break;
		}
	}
}

G3SkyMapPtr
FlatSkyMap::Clone(bool copy_data) const
{
	if (copy_data)
		return boost::make_shared<FlatSkyMap>(*this);
	else
		return boost::make_shared<FlatSkyMap>(proj_info,
		    coord_ref, is_weighted, units, pol_type);
}

double
FlatSkyMap::operator [] (int i) const
{
	if (dense_)
		return (*dense_)(i % xpix_, i / xpix_);
	if (sparse_)
		return (*sparse_)(i % xpix_, i / xpix_);
	return 0;
}

double &
FlatSkyMap::operator [] (int i)
{
	assert(i >= 0);
	assert(i < xpix_*ypix_);

	if (dense_)
		return (*dense_)(i % xpix_, i / xpix_);
	if (!sparse_)
		sparse_ = new SparseMapData(xpix_, ypix_);
	return (*sparse_)(i % xpix_, i / xpix_);
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

#define GETSET(name, type)                     \
	type FlatSkyMap::name() const          \
	{                                      \
		return proj_info.name();       \
	}                                      \
	void FlatSkyMap::set_##name(type name) \
	{                                      \
		proj_info.set_##name(name);    \
	}

GETSET(proj, MapProjection);
GETSET(alpha_center, double);
GETSET(delta_center, double);
GETSET(xres, double);
GETSET(yres, double);
GETSET(res, double);

std::vector<double> FlatSkyMap::angle_to_xy(double alpha, double delta) const {
	return proj_info.angle_to_xy(alpha, delta);
}

std::vector<double> FlatSkyMap::xy_to_angle(double x, double y) const {
	return proj_info.xy_to_angle(x, y, false);
}

size_t FlatSkyMap::angle_to_pixel(double alpha, double delta) const {
	return proj_info.angle_to_pixel(alpha, delta);
}

std::vector<double> FlatSkyMap::pixel_to_angle(size_t pixel) const {
	return proj_info.pixel_to_angle(pixel, false);
}

std::vector<double> FlatSkyMap::pixel_to_angle_wrap_ra(size_t pixel) const {
	return proj_info.pixel_to_angle(pixel, true);
}

void FlatSkyMap::get_rebin_angles(long pixel, size_t scale,
    std::vector<double> & alphas, std::vector<double> & deltas) const
{
	proj_info.get_rebin_angles(pixel, scale, alphas, deltas, false);
}

void FlatSkyMap::get_interp_pixels_weights(double alpha, double delta,
    std::vector<long> & pixels, std::vector<double> & weights) const
{
	proj_info.get_interp_pixels_weights(alpha, delta, pixels, weights);
}

G3SkyMapPtr FlatSkyMap::rebin(size_t scale) const
{
	if ((xpix_ % scale != 0) || (ypix_ % scale != 0)) {
		log_fatal("Map dimensions must be a multiple of rebinning scale");
	}

	if (scale <= 1)
		return Clone(true);

	FlatSkyProjection p(proj_info.rebin(scale));
	FlatSkyMapPtr out(new FlatSkyMap(p, coord_ref, is_weighted, units, pol_type));

#if 0
	if (IsAllocated()) {
		out->EnsureAllocated();

		for (size_t i = 0; i < xpix_; i++) {
			for (size_t j = 0; j < ypix_; j++) {
				(*out)[out->pixat(i / scale, j / scale)] += data_[pixat(i, j)];
			}
		}

		out /= (scale * scale);
	}
#endif
	return out;
}


#define FLAT_SKY_MAP_DOCSTR \
        "FlatSkyMap is a G3SkyMap with the extra meta information about the" \
	" particular flat sky projection included.  In practice it behaves\n" \
	" (mostly) like a 2d numpy array.\n\n"				\
        "Meta Information Stored: \n\n" \
        "    x_len (int) x (ra/az)  dimension length \n" \
        "    y_len (int) y (dec/el) dimension  length \n" \
        "    alpha_center : (double)  Ra (or Az) of the center of the map \n" \
        "    delta_center : (double)  Dec (or El) of the center of the map \n" \
        "    res : (double) approximate resolution of the pixel\n" \
	"          (the actual shape of a pixel is projection dependent)\n"\
        "    proj : (MapProjection)  proj is a MapProjection enum that specifies\n"\
	"           the flat sky projection \n"	\
        "\n\n" \
        "The other meta information is inherited from G3SkyMap that lives in core. \n\n" \
        "For reasons (skymap __setitem__ has to handle both 1d and 2d \n"\
	" semantics) the FlatSkyMap has a slightly unintuitive way of \n"\
	" setting values when using a slicing operator.  Instead of being\n"\
	" able to  slice directly you need to cast it to be an array first: \n\n" \
	"    np.asarray(your_flat_sky_map)[:] = the_numpy_array_you_are_assigning\n\n\n"\
	"If you find that you need numpy functionality from a FlatSkyMap,\n"\
	" using np.asarray will convert it to a numpy array without copying the data.\n" \
	" any changes to the resulting numpy array will affect the data stored in the\n" \
	" FlatSkyMap."



G3_SPLIT_SERIALIZABLE_CODE(FlatSkyMap);

PYBINDINGS("coordinateutils")
{
	using namespace boost::python;

	// Can't use the normal FRAMEOBJECT code since this inherits
	// from an intermediate class. Expanded by hand here.
	class_<FlatSkyMap, bases<G3SkyMap>, FlatSkyMapPtr>(
	  "FlatSkyMap", FLAT_SKY_MAP_DOCSTR, boost::python::no_init)
	    .def(boost::python::init<const FlatSkyMap &>())
	    .def_pickle(g3frameobject_picklesuite<FlatSkyMap>())
	    .def(bp::init<int, int, double, bool, MapProjection, double,
	       double, MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType, double>(
	         (bp::arg("x_len"), bp::arg("y_len"), bp::arg("res"),
		  bp::args("is_weighted") = true,
		  bp::arg("proj") = MapProjection::ProjNone,
		  bp::arg("alpha_center") = 0, bp::arg("delta_center") = 0,
		  bp::arg("coord_ref") = MapCoordReference::Equatorial,
		  bp::arg("units") = G3Timestream::Tcmb,
		  bp::arg("pol_type") = G3SkyMap::None, bp::arg("x_res") = 0)))
	    .def(bp::init<boost::python::object, double, bool, MapProjection,
	       double, double, MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType, double>(
		  (bp::arg("obj"), bp::arg("res"),
	           bp::args("is_weighted") = true,
		   bp::arg("proj") = MapProjection::ProjNone,
		   bp::arg("alpha_center") = 0, bp::arg("delta_center") = 0,
		   bp::arg("coord_ref") = MapCoordReference::Equatorial,
		   bp::arg("units") = G3Timestream::Tcmb,
		   bp::arg("pol_type") = G3SkyMap::None,
		   bp::arg("x_res") = 0)))

	    .def(bp::init<const FlatSkyMap&>(bp::arg("flat_map")))
	    .def(bp::init<>())
	    .add_property("proj", &FlatSkyMap::proj, &FlatSkyMap::set_proj,
	      "Map projection (one of coordinateutils.MapProjection)")
	    .add_property("alpha_center", &FlatSkyMap::alpha_center,
	      &FlatSkyMap::set_alpha_center, "Horizontal axis center position")
	    .add_property("delta_center", &FlatSkyMap::delta_center,
	      &FlatSkyMap::set_delta_center, "Vertical axis center position")
	    .add_property("res", &FlatSkyMap::res, &FlatSkyMap::set_res,
	      "Map resolution in angular units for maps with square pixels")
	    .add_property("x_res", &FlatSkyMap::xres, &FlatSkyMap::set_xres,
	      "Resolution in X direction for maps with rectangular pixels")
	    .add_property("y_res", &FlatSkyMap::yres, &FlatSkyMap::set_yres,
	      "Resolution in Y direction for maps with rectangular pixels")

	    .def("xy_to_angle", &FlatSkyMap::xy_to_angle,
	      (bp::arg("x"), bp::arg("y")),
	       "Compute the sky coordinates of the input flat 2D coordinates")
	    .def("angle_to_xy", &FlatSkyMap::angle_to_xy,
	      (bp::arg("alpha"), bp::arg("delta")),
	       "Compute the flat 2D coordinates of the input sky coordinates")
	;
	register_pointer_conversions<FlatSkyMap>();

	implicitly_convertible<FlatSkyMapPtr, G3SkyMapPtr>();
	implicitly_convertible<FlatSkyMapConstPtr, G3SkyMapConstPtr>();
	implicitly_convertible<FlatSkyMapPtr, G3SkyMapConstPtr>();
	implicitly_convertible<FlatSkyMapConstPtr, G3SkyMapConstPtr>();
}

