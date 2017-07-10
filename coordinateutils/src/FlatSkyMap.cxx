#include <pybindings.h>
#include <serialization.h>

#include <coordinateutils/FlatSkyMap.h>

FlatSkyMap::FlatSkyMap(int x_len, int y_len, double res, bool is_weighted,
    MapProjection proj, double alpha_center, double delta_center,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, double x_res) :
      G3SkyMap(coord_ref, x_len, y_len, is_weighted, u, pol_type, false),
      proj(proj), alpha_center(alpha_center), delta_center(delta_center),
      res(res), x_res(x_res > 0 ? x_res : res)
{
}

// Needs to pass-through to G3SkyMap constructor, so can't be an
// out-of-class make_constructor() thing
FlatSkyMap::FlatSkyMap(boost::python::object v, double res,
    bool is_weighted, MapProjection proj,
    double alpha_center, double delta_center,
    MapCoordReference coord_ref, G3Timestream::TimestreamUnits u,
    G3SkyMap::MapPolType pol_type, double x_res) :
      G3SkyMap(v, coord_ref, is_weighted, u, pol_type),
      proj(proj), alpha_center(alpha_center), delta_center(delta_center),
      res(res), x_res(x_res > 0 ? x_res : res)
{
}

FlatSkyMap::FlatSkyMap() :
    G3SkyMap(MapCoordReference::Local, 0), proj(MapProjection::Proj0),
    alpha_center(0), delta_center(0), res(0), x_res(0)
{
}

FlatSkyMap::FlatSkyMap(const FlatSkyMap & fm) :
    G3SkyMap(fm), proj(fm.proj), alpha_center(fm.alpha_center),
    delta_center(fm.delta_center), res(fm.res), x_res(fm.x_res)
{
}

template <class A> void
FlatSkyMap::serialize(A &ar, const unsigned u)
{
	using namespace cereal;

	ar & make_nvp("G3SkyMap", base_class<G3SkyMap>(this));
	ar & make_nvp("proj",proj);
	ar & make_nvp("alpha_center",alpha_center);
	ar & make_nvp("delta_center",delta_center);
	ar & make_nvp("res",res);
	ar & make_nvp("x_res", x_res);
}

G3SkyMapPtr
FlatSkyMap::Clone(bool copy_data) const
{
	if (copy_data)
		return boost::make_shared<FlatSkyMap>(*this);
	else
		return boost::make_shared<FlatSkyMap>(xpix_,  ypix_, res,
		    is_weighted, proj, alpha_center, delta_center, coord_ref,
		    units, pol_type, x_res);
}

std::string
FlatSkyMap::Description() const
{
	std::ostringstream os;

	os.precision(1);

	os << xpix_ << " x " << ypix_ <<
	    " (" << xpix_*((x_res != 0) ? x_res : res)/G3Units::deg << " x "
	    " (" << ypix_*res/G3Units::deg << " deg) in ";

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
	os << " coordinates, ";

	switch (proj) {
	case ProjSansonFlamsteed:
		os << "Sanson-Flamsteed";
		break;
	case ProjCAR:
		os << "CAR";
		break;
	case ProjSIN:
		os << "Sin";
		break;
	case ProjStereographic:
		os << "stereographic";
		break;
	case ProjLambertAzimuthalEqualArea:
		os << "Lambert Azimuthal Equal Area";
		break;
	case ProjBICEP:
		os << "BICEP";
		break;
	default:
		os << "other (" << proj << ")";
	}

	switch (units) {
	case G3Timestream::Counts:
		os << " (Counts)";
		break;
	case G3Timestream::Amps:
		os << " (Amps)";
		break;
	case G3Timestream::Watts:
		os << " (Watts)";
		break;
	case G3Timestream::Ohms:
		os << " (Ohms)";
		break;
	case G3Timestream::Kcmb:
		os << " (Kcmb)";
		break;
	default:
		break;
	}

	return os.str();
}

FlatSkyMap &
FlatSkyMap::operator+=(const FlatSkyMap & rhs)
{
	EnsureAllocated();
	g3_assert(data_.size() == rhs.data_.size());
	for (size_t i = 0; i < data_.size(); i++)
		data_[i] += rhs.data_[i];
	return *this;
}

FlatSkyMap &
FlatSkyMap::operator+=(double rhs)
{
	EnsureAllocated();
	for (size_t i = 0; i < data_.size(); i++)
		data_[i] += rhs;
	return *this;
}

FlatSkyMap
FlatSkyMap::operator+(const FlatSkyMap & rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map += rhs;
	return new_map;
}

FlatSkyMap
FlatSkyMap::operator+(double rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map += rhs;
	return new_map;
}

FlatSkyMap &
FlatSkyMap::operator-=(const FlatSkyMap &rhs)
{
	EnsureAllocated();
	g3_assert(data_.size() == rhs.data_.size());
	for (size_t i = 0; i < data_.size(); i++)
		data_[i] -= rhs.data_[i];
	return *this;
}

FlatSkyMap &
FlatSkyMap::operator-=(double rhs)
{
	EnsureAllocated();
	for (size_t i = 0; i < data_.size(); i++)
	data_[i] -= rhs;
	return *this;
}

FlatSkyMap
FlatSkyMap::operator-(const FlatSkyMap & rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map -= rhs;
	return new_map;
}

FlatSkyMap
FlatSkyMap::operator-(double rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map -= rhs;
	return new_map;
}

FlatSkyMap &
FlatSkyMap::operator*=(const FlatSkyMap & rhs)
{
	EnsureAllocated();
	g3_assert(data_.size() == rhs.data_.size());
	for (size_t i = 0; i < data_.size(); i++)
		data_[i] *= rhs.data_[i];
	return *this;
}

FlatSkyMap &
FlatSkyMap::operator*=(double rhs)
{
	EnsureAllocated();
	for (size_t i = 0; i < data_.size(); i++)
		data_[i] *= rhs;
	return *this;
}

FlatSkyMap
FlatSkyMap::operator*(const FlatSkyMap & rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map *= rhs;
	return new_map;
}

FlatSkyMap
FlatSkyMap::operator*(double rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map *= rhs;
	return new_map;
}

FlatSkyMap &
FlatSkyMap::operator/=(const FlatSkyMap & rhs)
{
	EnsureAllocated();
	g3_assert(data_.size() == rhs.data_.size());
	for (size_t i = 0; i < data_.size(); i++)
		data_[i] /= rhs.data_[i];
	return *this;
}

FlatSkyMap &
FlatSkyMap::operator/=(double rhs)
{
	EnsureAllocated();
	for (size_t i = 0; i < data_.size(); i++)
		data_[i] /= rhs;
	return *this;
}

FlatSkyMap
FlatSkyMap::operator/(const FlatSkyMap & rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map /= rhs;
	return new_map;
}

FlatSkyMap
FlatSkyMap::operator/(double rhs)
{
	EnsureAllocated();
	FlatSkyMap new_map(*this);
	new_map /= rhs;
	return new_map;
}

#define FLAT_SKY_MAP_DOCSTR \
        "FlatSkyMap is a G3SkyMap with the extra meta information about the particular flat sky projection included\n\n" \
        "For reasons (skymap __setitem__ has to handle both 1d and 2d semantics) the FlatSkyMap has a slightly unintuitive way of setting values when using a slicing operator.  Instead of being able to  slice directly you need to cast it to be an array first:\n\n"\
	"    np.asarray(your_flat_sky_map)[:] = the_numpy_array_you_are_assigning\n\n\n"\
        "Meta Information Stored: \n\n" \
        "    x_len (int) x (ra/az)  dimension length \n" \
        "    y_len (int) y (dec/el) dimension  length \n" \
        "    alpha_center : (double)  Ra (or Az) of the center of the map \n" \
        "    delta_center : (double)  Dec (or El) of the center of the map \n" \
        "    res : (double) approximate resolution of the pixel (this is projection dependent) \n" \
        "    proj : (MapProjection)  proj is a MapProjection enum that specifies the projection \n" \
        "\n\n" \
        "The other meta information is inheritted from G3SkyMap that lives in core. \n" \

G3_SERIALIZABLE_CODE(FlatSkyMap);

PYBINDINGS("coordinateutils")
{
	using namespace boost::python;

	bp::enum_<MapProjection>("MapProjection")
	    .value("Proj0", Proj0)
	    .value("Proj1", Proj1)
	    .value("Proj2", Proj2)
	    .value("Proj3", Proj3)
	    .value("Proj4", Proj4)
	    .value("Proj5", Proj5)
	    .value("Proj6", Proj6)
	    .value("Proj7", Proj7)
	    .value("Proj8", Proj8)
	    .value("Proj9", Proj9)

	    .value("ProjSansonFlamsteed", ProjSansonFlamsteed)
	    .value("ProjCAR", ProjCAR)
	    .value("ProjSIN", ProjSIN)
	    .value("ProjStereographic", ProjStereographic)
	    .value("ProjLambertAzimuthalEqualArea",
	      ProjLambertAzimuthalEqualArea)
	    .value("ProjBICEP", ProjBICEP)
	    .value("ProjNone", ProjNone)
	;

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
		  bp::arg("units") = G3Timestream::Kcmb,
		  bp::arg("pol_type") = G3SkyMap::None, bp::arg("x_res") = 0)))
	    .def(bp::init<boost::python::object, double, bool, MapProjection,
	       double, double, MapCoordReference, G3Timestream::TimestreamUnits,
	       G3SkyMap::MapPolType, double>(
		  (bp::arg("obj"), bp::arg("res"),
	           bp::args("is_weighted") = true,
		   bp::arg("proj") = MapProjection::ProjNone,
		   bp::arg("alpha_center") = 0, bp::arg("delta_center") = 0,
		   bp::arg("coord_ref") = MapCoordReference::Equatorial,
		   bp::arg("units") = G3Timestream::Kcmb,
		   bp::arg("pol_type") = G3SkyMap::None,
		   bp::arg("x_res") = 0)))

	    .def(bp::init<const FlatSkyMap&>(bp::arg("flat_map")))
	    .def(bp::init<>())
	    .def_readwrite("proj", &FlatSkyMap::proj,
	      "Map projection (one of coordinateutils.MapProjection)")
	    .def_readwrite("alpha_center", &FlatSkyMap::alpha_center,
	      "Horizontal axis center position")
	    .def_readwrite("delta_center", &FlatSkyMap::delta_center,
	      "Vertical axis center position")
	    .def_readwrite("res", &FlatSkyMap::res, "Map resolution in "
	      "angular units for maps with square pixels")
	    .def_readwrite("x_res", &FlatSkyMap::x_res, "Resolution in X "
	      "direction for maps with rectangular pixels")
	    .def_readwrite("y_res", &FlatSkyMap::res, "Resolution in Y "
	      "direction for maps with rectangular pixels")

	    .def(bp::self + bp::self)
	    .def(bp::self * bp::self)
	    .def(bp::self - bp::self)
	    .def(bp::self / bp::self)
	    .def(bp::self + double())
	    .def(bp::self * double())
	    .def(bp::self - double())
	    .def(bp::self / double())
	    .def(bp::self += bp::self)
	    .def(bp::self *= bp::self)
	    .def(bp::self -= bp::self)
	    .def(bp::self /= bp::self)
	    .def(bp::self += double())
	    .def(bp::self *= double())
	    .def(bp::self -= double())
	    .def(bp::self /= double())
	    .def("Clone", &FlatSkyMap::Clone)
	    .def("__copy__", &FlatSkyMap::Clone)
	;
	register_pointer_conversions<FlatSkyMap>();

	implicitly_convertible<FlatSkyMapPtr, G3SkyMapPtr>();
	implicitly_convertible<FlatSkyMapConstPtr, G3SkyMapConstPtr>();
	implicitly_convertible<FlatSkyMapPtr, G3SkyMapConstPtr>();
	implicitly_convertible<FlatSkyMapConstPtr, G3SkyMapConstPtr>();
}

