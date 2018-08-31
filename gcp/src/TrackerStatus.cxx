#include <pybindings.h>
#include <serialization.h>
#include <gcp/TrackerStatus.h>

template <class A> void TrackerStatus::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("time", time);
	ar & make_nvp("az_pos", az_pos);
	ar & make_nvp("el_pos", el_pos);
	ar & make_nvp("az_rate", az_rate);
	ar & make_nvp("el_rate", el_rate);
	ar & make_nvp("az_command", az_command);
	ar & make_nvp("el_command", el_command);
	ar & make_nvp("az_rate_command", az_rate_command);
	ar & make_nvp("el_rate_command", el_rate_command);

	ar & make_nvp("state", state);
	ar & make_nvp("acu_seq", acu_seq);
	ar & make_nvp("in_control", in_control);
	ar & make_nvp("scan_flag", scan_flag);
}

TrackerStatus TrackerStatus::operator +(const TrackerStatus &a) const
{
	TrackerStatus t(*this);

	t += a;
	return t;
}

TrackerStatus &TrackerStatus::operator +=(const TrackerStatus &a)
{
	#define tracker_appendarr(x) x.insert(x.end(), a.x.begin(), a.x.end())

	tracker_appendarr(time);
	tracker_appendarr(az_pos);
	tracker_appendarr(el_pos);
	tracker_appendarr(az_rate);
	tracker_appendarr(el_rate);
	tracker_appendarr(az_command);
	tracker_appendarr(el_command);
	tracker_appendarr(az_rate_command);
	tracker_appendarr(el_rate_command);
	tracker_appendarr(state);
	tracker_appendarr(acu_seq);
	tracker_appendarr(in_control);
	tracker_appendarr(scan_flag);

	#undef tracker_appendarr

	return *this;
}

std::string TrackerStatus::Description() const
{
	std::ostringstream s;

	s << time.size() << " tracker samples";

	if (time.size() > 0)
		s << " from " << time[0] << " to " << time[time.size() - 1];

	return s.str();
}

G3_SERIALIZABLE_CODE(TrackerStatus);

PYBINDINGS("gcp") {
	using namespace boost::python;

	enum_<enum TrackerStatus::TrackerState>("TrackerState")
	    .value("LACKING", TrackerStatus::LACKING)
	    .value("TIME_ERROR", TrackerStatus::TIME_ERROR)
	    .value("UPDATING", TrackerStatus::UPDATING)
	    .value("HALTED", TrackerStatus::HALTED)
	    .value("SLEWING", TrackerStatus::SLEWING)
	    .value("TRACKING", TrackerStatus::TRACKING)
	    .value("TOO_LOW", TrackerStatus::TOO_LOW)
	    .value("TOO_HIGH", TrackerStatus::TOO_HIGH)
	;
	register_vector_of<enum TrackerStatus::TrackerState>("TrackerState");

	EXPORT_FRAMEOBJECT(TrackerStatus, init<>(), "GCP Tracker Status")
	    .def_readwrite("time", &TrackerStatus::time)
	    .def_readwrite("az_pos", &TrackerStatus::az_pos)
	    .def_readwrite("el_pos", &TrackerStatus::el_pos)
	    .def_readwrite("az_rate", &TrackerStatus::az_rate)
	    .def_readwrite("el_rate", &TrackerStatus::el_rate)
	    .def_readwrite("az_command", &TrackerStatus::az_command)
	    .def_readwrite("el_command", &TrackerStatus::el_command)
	    .def_readwrite("az_rate_command", &TrackerStatus::az_rate_command)
	    .def_readwrite("el_rate_command", &TrackerStatus::el_rate_command)
	    .def_readwrite("state", &TrackerStatus::state)
	    .def_readwrite("acu_seq", &TrackerStatus::acu_seq)
	    .def_readwrite("in_control", &TrackerStatus::in_control)
	    .def_readwrite("scan_flag", &TrackerStatus::scan_flag)
	    .def(self + self)
	    .def(self += self)
	;
}

