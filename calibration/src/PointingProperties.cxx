#include <pybindings.h>
#include <calibration/PointingProperties.h>
#include <container_pybindings.h>

template <class A> void PointingProperties::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("tiltLat", tiltLat);
	ar & make_nvp("tiltHA", tiltHA);
	ar & make_nvp("tiltMag", tiltMag);
	ar & make_nvp("tiltAngle", tiltAngle);
}

std::string PointingProperties::Description() const
{
	std::ostringstream s;
	s << "Pointing model properties";

	return s.str();
}

G3_SERIALIZABLE_CODE(PointingProperties);
G3_SERIALIZABLE_CODE(PointingPropertiesMap);

PYBINDINGS("calibration", scope) {
	register_frameobject<PointingProperties>(scope, "PointingProperties",
            "Pointing model parameters to be used for offline pointing corrections.")
	    .def(py::init<>())
	    .def_readwrite("tiltLat", &PointingProperties::tiltLat,
	       "Azimuth lateral tilt parameter.")
	    .def_readwrite("tiltHA", &PointingProperties::tiltHA,
	       "Azimuth hour angle tilt parameter.")
	    .def_readwrite("tiltMag", &PointingProperties::tiltMag,
	       "Magnitude of azimuth tilt.")
	    .def_readwrite("tiltAngle", &PointingProperties::tiltAngle,
	       "Orientation of azimuth tilt.")
	;

	register_g3map<PointingPropertiesMap>("PointingPropertiesMap",
	    "Container for pointing model parameters for offline pointing.");
}

