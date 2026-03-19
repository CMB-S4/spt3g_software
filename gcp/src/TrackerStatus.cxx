#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
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

	if (v > 1) {
		ar & make_nvp("lst", lst);
		ar & make_nvp("source_acquired", source_acquired);
		ar & make_nvp("source_acquired_threshold", source_acquired_threshold);
		ar & make_nvp("tracker_mode", tracker_mode);
		ar & make_nvp("tracker_lacking", tracker_lacking);
		ar & make_nvp("time_status", time_status);
		ar & make_nvp("schedule_name", schedule_name);

		in_control_int = std::vector<int>(in_control.begin(), in_control.end());
	}
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

	tracker_appendarr(in_control_int);
	tracker_appendarr(lst);
	tracker_appendarr(source_acquired);
	tracker_appendarr(source_acquired_threshold);
	tracker_appendarr(tracker_mode);
	tracker_appendarr(tracker_lacking);
	tracker_appendarr(time_status);
	tracker_appendarr(schedule_name);

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

PYBINDINGS("gcp", scope) {
	register_enum<TrackerStatus::TrackerState>(scope, "TrackerState")
	    .value("LACKING", TrackerStatus::LACKING)
	    .value("TIME_ERROR", TrackerStatus::TIME_ERROR)
	    .value("UPDATING", TrackerStatus::UPDATING)
	    .value("HALTED", TrackerStatus::HALTED)
	    .value("SLEWING", TrackerStatus::SLEWING)
	    .value("TRACKING", TrackerStatus::TRACKING)
	    .value("TOO_LOW", TrackerStatus::TOO_LOW)
	    .value("TOO_HIGH", TrackerStatus::TOO_HIGH)
	;
	register_vector_of<TrackerStatus::TrackerState>(scope, "TrackerState");

	register_frameobject<TrackerStatus>(scope, "TrackerStatus", "GCP Tracker Status")
	    .def(py::init<>())
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
	    .def_readwrite("in_control_int", &TrackerStatus::in_control_int)
	    .def_readwrite("lst", &TrackerStatus::lst)
	    .def_readwrite("source_acquired", &TrackerStatus::source_acquired)
	    .def_readwrite("source_acquired_threshold", &TrackerStatus::source_acquired_threshold)
	    .def_readwrite("tracker_mode", &TrackerStatus::tracker_mode)
	    .def_readwrite("tracker_lacking", &TrackerStatus::tracker_lacking)
	    .def_readwrite("time_status", &TrackerStatus::time_status)
	    .def_readwrite("schedule_name", &TrackerStatus::schedule_name)
	    .def(py::self + py::self)
	    .def(py::self += py::self)
	;
}

