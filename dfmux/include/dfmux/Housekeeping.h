#ifndef _DFMUX_HOUSEKEEPING_H
#define _DFMUX_HOUSEKEEPING_H

#include <G3Frame.h>
#include <G3TimeStamp.h>
#include <G3Map.h>

#include <stdint.h>

/*
 * DfMux housekeeping data is organized into a tree:
 * HkBoardInfo -> HkMezzanineInfo ->  HkModuleInfo -> HkChannelInfo
 * (board level)  (mezz / SQCB)       (DAC/ADC/SQUID) (channel)
 *
 * Housekeeping data as broadly defined is both board configuration
 * and status (rails, temperatures, voltages). Please see IceBoard
 * documentation for details on the meaning of specific parameters.
 * Note that, unlike the wiring map, channels, modules, and mezzanines
 * are *1-indexed* when they occur in housekeeping data.
 */

class HkChannelInfo : public G3FrameObject
{
public:
	HkChannelInfo() : G3FrameObject(), channel_number(-1), carrier_amplitude(NAN),
        carrier_frequency(NAN), dan_accumulator_enable(false), dan_feedback_enable(false),
        dan_streaming_enable(false), dan_gain(NAN), demod_frequency(NAN), nuller_amplitude(NAN),
        dan_railed(false), rlatched(NAN), rnormal(NAN), rfrac_achieved(NAN), loopgain(NAN) {}

	int32_t channel_number; // 1-indexed
	double carrier_amplitude;
	double carrier_frequency;
	bool dan_accumulator_enable;
	bool dan_feedback_enable;
	bool dan_streaming_enable;
	double dan_gain;
	double demod_frequency;
	double nuller_amplitude;

	bool dan_railed;

	std::string state;
	double rlatched;
	double rnormal;
	double rfrac_achieved;
	double loopgain;

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
};

G3_POINTERS(HkChannelInfo);
G3_SERIALIZABLE(HkChannelInfo, 5);

class HkModuleInfo : public G3FrameObject
{
public:
	HkModuleInfo() : G3FrameObject(), module_number(-1), carrier_gain(-1),
        nuller_gain(-1), demod_gain(-1), carrier_railed(false), nuller_railed(false),
        demod_railed(false), squid_flux_bias(NAN), squid_current_bias(NAN),
        squid_stage1_offset(NAN), squid_p2p(NAN), squid_transimpedance(NAN) {}

	int32_t module_number; // 1-indexed

	int32_t carrier_gain;
	int32_t nuller_gain;
	int32_t demod_gain;

	bool carrier_railed;
	bool nuller_railed;
	bool demod_railed;

	double squid_flux_bias;
	double squid_current_bias;
	double squid_stage1_offset;

	double squid_p2p;
	double squid_transimpedance;
	std::string squid_state;

	std::string squid_feedback;
	std::string routing_type;

	std::map<int32_t, HkChannelInfo> channels; // 1-indexed

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
};

G3_POINTERS(HkModuleInfo);
G3_SERIALIZABLE(HkModuleInfo, 2);

class HkMezzanineInfo : public G3FrameObject
{
public:
	HkMezzanineInfo() : G3FrameObject(), power(false), present(false),
       temperature(NAN), squid_controller_temperature(NAN), squid_heater(NAN),
       squid_controller_power(false) {}

	bool power;
	bool present;

	std::string serial;
	std::string part_number;
	std::string revision;

	std::map<std::string, double> currents;
	std::map<std::string, double> voltages;

	std::map<int32_t, HkModuleInfo> modules; // 1-indexed

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;

	double temperature;
	double squid_controller_temperature;
	double squid_heater;
	bool squid_controller_power;
};

G3_POINTERS(HkMezzanineInfo);
G3_SERIALIZABLE(HkMezzanineInfo, 2);

class HkBoardInfo : public G3FrameObject
{
public:
	HkBoardInfo() : G3FrameObject(), fir_stage(-1), is128x(false) {}

	G3Time timestamp;

	std::string timestamp_port;
	std::string serial;
	int32_t fir_stage;
	bool is128x;

	std::map<std::string, double> currents;
	std::map<std::string, double> voltages;
	std::map<std::string, double> temperatures;

	std::map<int32_t, HkMezzanineInfo> mezz; // 1-indexed

	template <class A> void serialize(A &ar, unsigned v);
	std::string Description() const;
};

G3_POINTERS(HkBoardInfo);
G3_SERIALIZABLE(HkBoardInfo, 2);

G3MAP_OF(int32_t, HkBoardInfo, DfMuxHousekeepingMap);

#endif
