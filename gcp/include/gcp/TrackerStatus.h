#ifndef TRACKERSTATUS_H
#define TRACKERSTATUS_H

#include <G3Frame.h>
#include <G3Timestream.h>
#include <string>

/*
 * Class containg Tracker information from GCP
 */
class TrackerStatus : public G3FrameObject
{
public:
	TrackerStatus() : G3FrameObject() {}
	std::string Description() const;

	// Concatenate two sets of arrays together
	TrackerStatus operator + (const TrackerStatus &) const;
	TrackerStatus &operator += (const TrackerStatus &);

	std::vector<G3Time> time;

	std::vector<double> az_pos, el_pos;
	std::vector<double> az_rate, el_rate;

	std::vector<double> az_command, el_command;
	std::vector<double> az_rate_command, el_rate_command;

	enum TrackerState {
		// See TrackingStatus.h in GCP for definitions
		LACKING,
		TIME_ERROR,
		UPDATING,
		HALTED,
		SLEWING,
		TRACKING,
		TOO_LOW,
		TOO_HIGH
	};

	std::vector<enum TrackerState> state;
	std::vector<int> acu_seq;
	std::vector<bool> in_control;
	std::vector<bool> scan_flag;

	// XXX identify which other parameters are meaningful and should be
	// in this container (e.g. tiltmeter status should be elsewhere)

	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(TrackerStatus);
G3_SERIALIZABLE(TrackerStatus, 1);

#endif

