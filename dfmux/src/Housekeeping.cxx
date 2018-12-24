#include <sstream>
#include <pybindings.h>
#include <dfmux/Housekeeping.h>

#include <container_pybindings.h>

std::string HkChannelInfo::Description() const
{
	std::ostringstream s;
	s << "Channel " << channel_number << ", " <<
	    carrier_frequency/G3Units::MHz << " MHz (tuning: " << state << ")";
	return s.str();
}

template <class A> void HkChannelInfo::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("channel_number", channel_number);
	ar & make_nvp("carrier_amplitude", carrier_amplitude);
	ar & make_nvp("carrier_frequency", carrier_frequency);
	ar & make_nvp("dan_accumulator_enable", dan_accumulator_enable);
	ar & make_nvp("dan_feedback_enable", dan_feedback_enable);
	ar & make_nvp("dan_streaming_enable", dan_streaming_enable);
	ar & make_nvp("dan_gain", dan_gain);
	ar & make_nvp("demod_frequency", demod_frequency);
	ar & make_nvp("nuller_amplitude", nuller_amplitude);
	ar & make_nvp("dan_railed", dan_railed);

	if (v >= 2) {
		ar & make_nvp("state", state);
		ar & make_nvp("rlatched", rlatched);
		ar & make_nvp("rnormal", rnormal);
		ar & make_nvp("rfrac_achieved", rfrac_achieved);
	}

	if (v == 3) {
		// Absorb a short-lived calibration constant that lived in
		// the HK data to address software limitations in Lyrebird.
		double res_conversion_factor;
		ar & make_nvp("res_conversion_factor", res_conversion_factor);
	}

	if (v > 4) {
		ar & make_nvp("loopgain", loopgain);
	}
}

std::string HkModuleInfo::Description() const
{
	std::ostringstream s;
	s << "Module " << module_number << " (SQUID: " << squid_state << ")";
	return s.str();
}

template <class A> void HkModuleInfo::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("module_number", module_number);
	ar & make_nvp("carrier_gain", carrier_gain);
	ar & make_nvp("nuller_gain", nuller_gain);
	ar & make_nvp("demod_gain", demod_gain);
	ar & make_nvp("carrier_railed", carrier_railed);
	ar & make_nvp("nuller_railed", nuller_railed);
	ar & make_nvp("demod_railed", demod_railed);
	ar & make_nvp("squid_flux_bias", squid_flux_bias);
	ar & make_nvp("squid_current_bias", squid_current_bias);
	ar & make_nvp("squid_stage1_offset", squid_stage1_offset);
	ar & make_nvp("squid_feedback", squid_feedback);
	ar & make_nvp("routing_type", routing_type);
	ar & make_nvp("channels", channels);

	if (v >= 2) {
		ar & make_nvp("squid_state", squid_state);
		ar & make_nvp("squid_p2p", squid_p2p);
		ar & make_nvp("squid_transimpedance", squid_transimpedance);
	}
}

std::string HkMezzanineInfo::Description() const
{
	std::ostringstream s;
	s << "Mezzanine serial " << serial << " (" << part_number <<
	    ") powered " << (power ? "on" : "off") << " and " <<
	    (present ? "" : "not ") << "present";
	return s.str();
}

template <class A> void HkMezzanineInfo::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("power", power);
	ar & make_nvp("present", present);
	ar & make_nvp("serial", serial);
	ar & make_nvp("part_number", part_number);
	ar & make_nvp("revision", revision);
	ar & make_nvp("currents", currents);
	ar & make_nvp("voltages", voltages);
	ar & make_nvp("modules", modules);

	if (v >= 2) {
		ar & make_nvp("temperature", temperature);
		ar & make_nvp("squid_controller_temperature", squid_controller_temperature);
		ar & make_nvp("squid_controller_power", squid_controller_power);
		ar & make_nvp("squid_heater", squid_heater);
	}
}

std::string HkBoardInfo::Description() const
{
	std::ostringstream s;
	s << "Board serial " << serial << ", set to FIR " << fir_stage <<
	    ", at " << timestamp;
	return s.str();
}

template <class A> void HkBoardInfo::serialize(A &ar, unsigned v)
{
	using namespace cereal;

	G3_CHECK_VERSION(v);

	ar & make_nvp("G3FrameObject", base_class<G3FrameObject>(this));
	ar & make_nvp("timestamp", timestamp);
	ar & make_nvp("timestamp_port", timestamp_port);
	ar & make_nvp("serial", serial);
	ar & make_nvp("fir_stage", fir_stage);

	ar & make_nvp("currents", currents);
	ar & make_nvp("voltages", voltages);
	ar & make_nvp("temperatures", temperatures);
	ar & make_nvp("mezz", mezz);

	if (v >= 2) {
		ar & make_nvp("is128x", is128x);
	}
}

G3_SERIALIZABLE_CODE(HkChannelInfo);
G3_SERIALIZABLE_CODE(HkModuleInfo);
G3_SERIALIZABLE_CODE(HkMezzanineInfo);
G3_SERIALIZABLE_CODE(HkBoardInfo);
G3_SERIALIZABLE_CODE(DfMuxHousekeepingMap);

PYBINDINGS("dfmux") {
	EXPORT_FRAMEOBJECT(HkChannelInfo, init<>(), "Mux channel status "
	  "(configuration and sensors). Usually a part of an HkModuleInfo.")
	    .def_readwrite("channel_number", &HkChannelInfo::channel_number,
               "1-indexed channel number.")
	    .def_readwrite("carrier_amplitude",
	       &HkChannelInfo::carrier_amplitude,
               "Carrier amplitude in normalized units (0-1).")
	    .def_readwrite("carrier_frequency",
	       &HkChannelInfo::carrier_frequency,
	       "Carrier frequency in standard frequency units")
	    .def_readwrite("dan_accumulator_enable",
	       &HkChannelInfo::dan_accumulator_enable,
	       "True if DAN accumulator enabled")
	    .def_readwrite("dan_feedback_enable",
	       &HkChannelInfo::dan_feedback_enable,
	       "True if DAN control of the nuller is enabled")
	    .def_readwrite("dan_streaming_enable",
	       &HkChannelInfo::dan_streaming_enable,
	       "True if streamer packets are from DAN rather than demod")
	    .def_readwrite("dan_gain", &HkChannelInfo::dan_gain,
	       "DAN gain in board units")
	    .def_readwrite("demod_frequency", &HkChannelInfo::demod_frequency,
	       "Frequency of the demodulator in standard frequency units.")
	    .def_readwrite("nuller_amplitude", &HkChannelInfo::nuller_amplitude,
	       "Nuller amplitude in normalized units (0-1)")
	    .def_readwrite("dan_railed", &HkChannelInfo::dan_railed,
	       "True if DAN has railed.")
	    .def_readwrite("state", &HkChannelInfo::state,
	       "String code (\"latched\", \"tuned\" etc.) describing the state "
	       "of this detector stored by the control software")
	    .def_readwrite("rlatched", &HkChannelInfo::rlatched,
	       "Resistance of the detector when latched in standard impedance "
	       "units as stored by the control software tuning script.")
	    .def_readwrite("rnormal", &HkChannelInfo::rnormal,
	       "Resistance of the detector when normal in standard impedance "
	       "units as stored by the control software tuning script.")
	    .def_readwrite("rfrac_achieved", &HkChannelInfo::rfrac_achieved,
	       "Achieved resistance of the detector when tuned as a fraction "
	       "of rnormal, as stored by the control software tuning script.")
	    .def_readwrite("loopgain", &HkChannelInfo::loopgain,
	       "Measured loopgain of the detector as stored by the "
	       "control software tuning script.")
	;
	register_map<std::map<int, HkChannelInfo> >("HkChannelInfoMap",
	    "Mapping of channel number (1-indexed) to channel status "
	    "information");

	EXPORT_FRAMEOBJECT(HkModuleInfo, init<>(), "Mux module status")
	    .def_readwrite("module_number", &HkModuleInfo::module_number,
	       "1-indexed module number on this mezzanine")
	    .def_readwrite("carrier_gain", &HkModuleInfo::carrier_gain,
	       "Carrier gain code, in board-specific units")
	    .def_readwrite("nuller_gain", &HkModuleInfo::nuller_gain,
	       "Nuller gain code, in board-specific units")
	    .def_readwrite("demod_gain", &HkModuleInfo::demod_gain,
	       "Demod gain code, in board-specific units")
	    .def_readwrite("carrier_railed", &HkModuleInfo::carrier_railed,
	       "True if carrier has recently hit a DAC rail")
	    .def_readwrite("nuller_railed", &HkModuleInfo::nuller_railed,
	       "True if nuller has recently hit a DAC rail")
	    .def_readwrite("demod_railed", &HkModuleInfo::demod_railed,
	       "True if demod input has recently hit one of the ADC rails")
	    .def_readwrite("squid_flux_bias", &HkModuleInfo::squid_flux_bias,
	       "Flux bias, in board-specific units (XXX)")
	    .def_readwrite("squid_current_bias", &HkModuleInfo::squid_current_bias,
	       "SQUID current bias, in board-specific units (XXX)")
	    .def_readwrite("squid_stage1_offset", &HkModuleInfo::squid_stage1_offset,
	       "Offset voltage applied to SQUID output before first-stage amp")
	    .def_readwrite("squid_p2p", &HkModuleInfo::squid_p2p,
	       "Peak-to-peak voltage difference in V-Phi curve stored by "
	       "tuning script")
	    .def_readwrite("squid_transimpedance", &HkModuleInfo::squid_transimpedance,
	       "Measured SQUID transimpedance, in resistance units")
	    .def_readwrite("squid_state", &HkModuleInfo::squid_state,
	       "Descriptive string (e.g. \"Tuned\") stored by control system "
	       "tuning script to describe SQUID state")
	    .def_readwrite("squid_feedback", &HkModuleInfo::squid_feedback,
	       "SQUID feedback mechanism employed")
	    .def_readwrite("routing_type", &HkModuleInfo::routing_type,
	       "Whether DAC are routed directly to ADCs or to the cryostat")
	    .def_readwrite("channels", &HkModuleInfo::channels,
	       "Mapping from 1-indexed channel numbers to channel housekeeping "
	       "data")
	;
	register_map<std::map<int, HkModuleInfo> >("HkModuleInfoMap",
	    "Mapping from 1-indexed module numbers to module-specific "
	    "housekeeping data");

	EXPORT_FRAMEOBJECT(HkMezzanineInfo, init<>(), "Mux mezzanine status")
	    .def_readwrite("power", &HkMezzanineInfo::power, "True if on")
	    .def_readwrite("present", &HkMezzanineInfo::present, "True if exists")
	    .def_readwrite("serial", &HkMezzanineInfo::serial,
	        "Mezzanine serial number")
	    .def_readwrite("part_number", &HkMezzanineInfo::part_number,
	        "Mezzanine part ID (usually \"MGMEZZ04\")")
	    .def_readwrite("revision", &HkMezzanineInfo::revision,
	        "Mezzanine revision number")
	    .def_readwrite("currents", &HkMezzanineInfo::currents,
	        "Dictionary of measured currents on mezzanine")
	    .def_readwrite("voltages", &HkMezzanineInfo::voltages,
	        "Dictionary of measured voltages on mezzanine")
	    .def_readwrite("modules", &HkMezzanineInfo::modules,
	        "1-indexed list of housekeeping data from readout modules on "
	        "this mezzanine")
	    .def_readwrite("temperature", &HkMezzanineInfo::temperature,
	        "Mezzanine temperature (C)")
	    .def_readwrite("squid_controller_temperature",
	        &HkMezzanineInfo::squid_controller_temperature,
	        "Measured temperature of SQUID controller board (C)")
	    .def_readwrite("squid_controller_power",
	        &HkMezzanineInfo::squid_controller_power,
	        "True if SQUID controller board powered up")
	    .def_readwrite("squid_heater", &HkMezzanineInfo::squid_heater,
	        "Power level of SQUID header control")
	;
	register_map<std::map<int, HkMezzanineInfo> >("HkMezzanineInfoMap",
	    "1-indexed mapping of mezzanine ID to mezzanine-specific "
	    "housekeeping data");

	EXPORT_FRAMEOBJECT(HkBoardInfo, init<>(), "Mux board status. Includes "
	  "both configuration and sensor readings for board generic quantities "
	  "and a list of quantities for the mezzanines.")
	    .def_readwrite("timestamp", &HkBoardInfo::timestamp,
	       "Time at which housekeeping data collected")
	    .def_readwrite("timestamp_port", &HkBoardInfo::timestamp_port,
	       "Source of timestamps on board")
	    .def_readwrite("serial", &HkBoardInfo::serial,
	       "Board serial number")
	    .def_readwrite("fir_stage", &HkBoardInfo::fir_stage,
	       "Sample rate encoded as the \"FIR Stage\". Smaller numbers are "
	       "faster and grow by factors of two with each decrement")
	    .def_readwrite("is128x", &HkBoardInfo::is128x,
		"Boolean for whether 128x firmware is running")
	    .def_readwrite("currents", &HkBoardInfo::currents,
	       "Dictionary of data from on-board current sensors")
	    .def_readwrite("voltages", &HkBoardInfo::voltages,
	       "Dictionary of data from on-board voltage sensors")
	    .def_readwrite("temperatures", &HkBoardInfo::temperatures,
	       "Dictionary of data from on-board temperature sensors (C)")
	    .def_readwrite("mezz", &HkBoardInfo::mezz,
	       "1-indexed mapping from mezzanine ID to mezzanine-specific data")
	;

	register_g3map<DfMuxHousekeepingMap>("DfMuxHousekeepingMap",
	   "Container structure for housekeeping data from all DfMux boards, "
	   "indexed by board serial number.");
}
