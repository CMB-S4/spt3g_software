#ifndef ACUSTATUS_H
#define ACUSTATUS_H

#include <G3Frame.h>
#include <G3Vector.h>
#include <G3TimeStamp.h>
#include <string>

/*
 * Class containg ACU information from GCP
 */
class ACUStatus : public G3FrameObject
{
public:
	ACUStatus() : az_pos(NAN), el_pos(NAN), az_rate(NAN), el_rate(NAN)
	    {}
	std::string Description() const;

	G3Time time;

	double az_pos; /* Focal plane position in angular units */
	double el_pos;

	double az_rate; /* Angular units per second */
	double el_rate;

	int px_checksum_error_count;
	int px_resync_count;
	int px_resync_timeout_count;
	int px_timeout_count;
	int restart_count;
	bool px_resyncing;

	enum ACUState {
		IDLE,
		TRACKING,
		WAIT_RESTART,
		RESYNC
	};

	enum ACUState state;
	int acu_status; // XXX Add bitmask definitions

	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(ACUStatus);
G3_SERIALIZABLE(ACUStatus, 2);

G3VECTOR_OF(ACUStatus, ACUStatusVector);

#endif

