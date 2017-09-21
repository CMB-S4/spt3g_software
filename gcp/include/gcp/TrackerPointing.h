#ifndef TRACKERPOINTING_H
#define TRACKERPOINTING_H

#include <G3Frame.h>
#include <G3Timestream.h>
#include <string>

/*
 * Class containing Tracker pointing information from GCP
 */
class TrackerPointing : public G3FrameObject
{
public:
	TrackerPointing() : G3FrameObject() {}
	std::string Description() const;

	// Concatenate two sets of arrays together
	TrackerPointing operator + (const TrackerPointing &) const;
	TrackerPointing &operator += (const TrackerPointing &);

	//Define slow register int vectors.
	std::vector<G3Time> time;
	std::vector<int32_t> features;

	//Define fast-register double vectors.
	std::vector<double> horiz_mount_x, horiz_mount_y;
	std::vector<double> horiz_off_x, horiz_off_y;
	std::vector<double> linsens_avg_l1, linsens_avg_l2;
	std::vector<double> linsens_avg_r1, linsens_avg_r2;
	std::vector<double> scu_temp;
	std::vector<double> telescope_temp, telescope_pressure;
	std::vector<double> encoder_off_x, encoder_off_y;

	//Define fast-register int vectors.
	std::vector<double> tilts_x, tilts_y;

	//Define stored refraction corrections
	std::vector<double> refraction;

	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(TrackerPointing);
G3_SERIALIZABLE(TrackerPointing, 2);

#endif

