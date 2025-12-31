#include <pybindings.h>
#include <serialization.h>
#include <gcp/TrackerPointing.h>

template <class A> void TrackerPointing::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	// These correspond to GCP slow registers (sampled once per frame).
	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("time", time);
	ar & make_nvp("scu_temp", scu_temp);
	ar & make_nvp("features", features);
	
	//Breaking horiz and encoder into its constituent components.
	ar & make_nvp("encoder_off_x", encoder_off_x);
	ar & make_nvp("encoder_off_y", encoder_off_y);
	ar & make_nvp("horiz_mount_x", horiz_mount_x);
	ar & make_nvp("horiz_mount_y", horiz_mount_y);
	ar & make_nvp("horiz_off_x", horiz_off_x);
	ar & make_nvp("horiz_off_y", horiz_off_y);


	if (v < 2) {
		//Old versions recorded online pointing model parameters here
		std::vector<double> junk;
		ar & make_nvp("tilts", junk);
		ar & make_nvp("flexure", junk);
		ar & make_nvp("fixedCollimation", junk);
	}

	//Online refraction corrections
	ar & make_nvp("refraction", refraction);

	//Breaking tilts into its constituent components.
	ar & make_nvp("tilts_x", tilts_x);
	ar & make_nvp("tilts_y", tilts_y);

	//Record each linear sensor individually
	ar & make_nvp("linsens_avg_l1", linsens_avg_l1);
	ar & make_nvp("linsens_avg_l2", linsens_avg_l2);
	ar & make_nvp("linsens_avg_r1", linsens_avg_r1);
	ar & make_nvp("linsens_avg_r2", linsens_avg_r2);
	
	//Record telescope temperature and pressure for refraction
	ar & make_nvp("telescope_temp", telescope_temp);
	ar & make_nvp("telescope_pressure", telescope_pressure);

	if (v > 2) {
		ar & make_nvp("low_limit_az", low_limit_az);
		ar & make_nvp("high_limit_az", high_limit_az);
		ar & make_nvp("low_limit_el", low_limit_el);
		ar & make_nvp("high_limit_el", high_limit_el);
		ar & make_nvp("scan_off_x", scan_off_x);
		ar & make_nvp("scan_off_y", scan_off_y);
		ar & make_nvp("sky_off_x", sky_off_x);
		ar & make_nvp("sky_off_y", sky_off_y);
		ar & make_nvp("equat_off_x", equat_off_x);
		ar & make_nvp("equat_off_y", equat_off_y);
		ar & make_nvp("equat_geoc_ra", equat_geoc_ra);
		ar & make_nvp("equat_geoc_dec", equat_geoc_dec);
		ar & make_nvp("horiz_topo_az", horiz_topo_az);
		ar & make_nvp("horiz_topo_el", horiz_topo_el);
		ar & make_nvp("error_az", error_az);
		ar & make_nvp("error_el", error_el);
	}
}

TrackerPointing TrackerPointing::operator +(const TrackerPointing &a) const
{
	TrackerPointing t(*this);

	t += a;
	return t;
}

TrackerPointing &TrackerPointing::operator +=(const TrackerPointing &a)
{
	#define tracker_appendarr(x) x.insert(x.end(), a.x.begin(), a.x.end())

        tracker_appendarr(time);
	tracker_appendarr(scu_temp);
	tracker_appendarr(features);
	tracker_appendarr(encoder_off_x);
	tracker_appendarr(encoder_off_y);
	tracker_appendarr(horiz_mount_x);
	tracker_appendarr(horiz_mount_y);
	tracker_appendarr(horiz_off_x);
	tracker_appendarr(horiz_off_y);
	tracker_appendarr(tilts_x);
	tracker_appendarr(tilts_y);
	tracker_appendarr(linsens_avg_l1);
	tracker_appendarr(linsens_avg_l2);
	tracker_appendarr(linsens_avg_r1);
	tracker_appendarr(linsens_avg_r2);
	tracker_appendarr(telescope_temp);
	tracker_appendarr(telescope_pressure);
	tracker_appendarr(refraction);

	tracker_appendarr(low_limit_az);
	tracker_appendarr(high_limit_az);
	tracker_appendarr(low_limit_el);
	tracker_appendarr(high_limit_el);
	tracker_appendarr(scan_off_x);
	tracker_appendarr(scan_off_y);
	tracker_appendarr(sky_off_x);
	tracker_appendarr(sky_off_y);
	tracker_appendarr(equat_off_x);
	tracker_appendarr(equat_off_y);
	tracker_appendarr(equat_geoc_ra);
	tracker_appendarr(equat_geoc_dec);
	tracker_appendarr(horiz_topo_az);
	tracker_appendarr(horiz_topo_el);
	tracker_appendarr(error_az);
	tracker_appendarr(error_el);

	#undef tracker_appendarr

	return *this;
}

std::string TrackerPointing::Description() const
{
	std::ostringstream s;

	s << time.size() << " tracker pointing samples";

	if (time.size() > 0)
		s << " from " << time[0] << " to " << time[time.size() - 1];

	return s.str();
}

G3_SERIALIZABLE_CODE(TrackerPointing);

PYBINDINGS("gcp", scope) {
	register_frameobject<TrackerPointing>(scope, "TrackerPointing", "GCP Tracker Pointing")
	    .def(py::init<>())
	    .def_readwrite("time", &TrackerPointing::time)
	    .def_readwrite("scu_temp", &TrackerPointing::scu_temp)
	    .def_readwrite("features", &TrackerPointing::features)
	    .def_readwrite("encoder_off_x", &TrackerPointing::encoder_off_x)
	    .def_readwrite("encoder_off_y", &TrackerPointing::encoder_off_y)
	    .def_readwrite("horiz_mount_x", &TrackerPointing::horiz_mount_x)
	    .def_readwrite("horiz_mount_y", &TrackerPointing::horiz_mount_y)
	    .def_readwrite("horiz_off_x", &TrackerPointing::horiz_off_x)
	    .def_readwrite("horiz_off_y", &TrackerPointing::horiz_off_y)
	    .def_readwrite("tilts_x", &TrackerPointing::tilts_x)
	    .def_readwrite("tilts_y", &TrackerPointing::tilts_y)
	    .def_readwrite("linsens_avg_l1", &TrackerPointing::linsens_avg_l1)
	    .def_readwrite("linsens_avg_l2", &TrackerPointing::linsens_avg_l2)
	    .def_readwrite("linsens_avg_r1", &TrackerPointing::linsens_avg_r1)
	    .def_readwrite("linsens_avg_r2", &TrackerPointing::linsens_avg_r2)
	    .def_readwrite("telescope_temp", &TrackerPointing::telescope_temp)
	    .def_readwrite("telescope_pressure", &TrackerPointing::telescope_pressure)
	    .def_readwrite("refraction", &TrackerPointing::refraction)
	    .def_readwrite("low_limit_az", &TrackerPointing::low_limit_az)
	    .def_readwrite("high_limit_az", &TrackerPointing::high_limit_az)
	    .def_readwrite("low_limit_el", &TrackerPointing::low_limit_el)
	    .def_readwrite("high_limit_el", &TrackerPointing::high_limit_el)
	    .def_readwrite("scan_off_x", &TrackerPointing::scan_off_x)
	    .def_readwrite("scan_off_y", &TrackerPointing::scan_off_y)
	    .def_readwrite("sky_off_x", &TrackerPointing::sky_off_x)
	    .def_readwrite("sky_off_y", &TrackerPointing::sky_off_y)
	    .def_readwrite("equat_off_x", &TrackerPointing::equat_off_x)
	    .def_readwrite("equat_off_y", &TrackerPointing::equat_off_y)
	    .def_readwrite("equat_geoc_ra", &TrackerPointing::equat_geoc_ra)
	    .def_readwrite("equat_geoc_dec", &TrackerPointing::equat_geoc_dec)
	    .def_readwrite("horiz_topo_az", &TrackerPointing::horiz_topo_az)
	    .def_readwrite("horiz_topo_el", &TrackerPointing::horiz_topo_el)
	    .def_readwrite("error_az", &TrackerPointing::error_az)
	    .def_readwrite("error_el", &TrackerPointing::error_el)
	    .def(py::self + py::self)
	    .def(py::self += py::self)
	;
}

