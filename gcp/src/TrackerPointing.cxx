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

PYBINDINGS("gcp") {
	using namespace boost::python;

	EXPORT_FRAMEOBJECT(TrackerPointing, init<>(), "GCP Tracker Pointing")
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
	    .def(self + self)
	    .def(self += self)
	;
}

