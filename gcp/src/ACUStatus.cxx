#include <pybindings.h>
#include <serialization.h>
#include <container_pybindings.h>
#include <gcp/ACUStatus.h>

template <class A> void ACUStatus::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("time", time);
	ar & make_nvp("az_pos", az_pos);
	ar & make_nvp("el_pos", el_pos);
	ar & make_nvp("az_rate", az_rate);
	ar & make_nvp("el_rate", el_rate);

	if (v < 2) {
		double az_err = 0, el_err = 0;
		ar & make_nvp("az_err", az_err);
		ar & make_nvp("el_err", el_err);
	}

	ar & make_nvp("px_checksum_error_count", px_checksum_error_count);
	ar & make_nvp("px_resync_count", px_resync_count);
	ar & make_nvp("px_resync_timeout_count", px_resync_timeout_count);
	ar & make_nvp("px_timeout_count", px_timeout_count);
	ar & make_nvp("restart_count", restart_count);
	ar & make_nvp("px_resyncing", px_resyncing);
	ar & make_nvp("state", state);
	ar & make_nvp("acu_status", acu_status);
}

std::string ACUStatus::Description() const
{
	std::ostringstream s;
	std::string state_string;

	switch (state) {
	case IDLE:
		state_string = "IDLE";
		break;
	case TRACKING:
		state_string = "TRACKING";
		break;
	case WAIT_RESTART:
		state_string = "WAIT RESTART";
		break;
	case RESYNC:
		state_string = "RESYNC";
		break;
	default:
		state_string = "Unknown ACU State";
	}

	s << "Az " << az_pos/G3Units::deg << " deg, el " << el_pos/G3Units::deg << " deg at " << time << ", " << state_string;

	return s.str();
}

G3_SERIALIZABLE_CODE(ACUStatus);
G3_SERIALIZABLE_CODE(ACUStatusVector);

static bool operator == (const ACUStatus &a, const ACUStatus &b) {
	throw std::runtime_error("ACU statuses cannot be compared");
	return false;
}

PYBINDINGS("gcp") {
	boost::python::enum_<enum ACUStatus::ACUState>("ACUState")
	    .value("IDLE", ACUStatus::IDLE)
	    .value("TRACKING", ACUStatus::TRACKING)
	    .value("WAIT_RESTART", ACUStatus::WAIT_RESTART)
	    .value("RESYNC", ACUStatus::RESYNC)
	;

	EXPORT_FRAMEOBJECT(ACUStatus, init<>(), "ACU Status information, as "
	  "reported by the ACU")
	    .def_readwrite("time", &ACUStatus::time)
	    .def_readwrite("az_pos", &ACUStatus::az_pos)
	    .def_readwrite("el_pos", &ACUStatus::el_pos)
	    .def_readwrite("az_rate", &ACUStatus::az_rate)
	    .def_readwrite("el_rate", &ACUStatus::el_rate)
	    .def_readwrite("px_checksum_error_count", &ACUStatus::px_checksum_error_count)
	    .def_readwrite("px_resync_count", &ACUStatus::px_resync_count)
	    .def_readwrite("px_resync_timeout_count", &ACUStatus::px_resync_timeout_count)
	    .def_readwrite("px_timeout_count", &ACUStatus::px_timeout_count)
	    .def_readwrite("restart_count", &ACUStatus::restart_count)
	    .def_readwrite("px_resyncing", &ACUStatus::px_resyncing)
	    .def_readwrite("state", &ACUStatus::state)
	    .def_readwrite("acu_status", &ACUStatus::acu_status)
	;

	register_vector_of<ACUStatus>("_ACUStatusVectorBase");
	register_g3vector<ACUStatus>("ACUStatusVector", "Array of ACUStatus "
	   "objects, usually time-ordered");
}

