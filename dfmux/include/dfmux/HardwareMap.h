#ifndef DFMUXHARDWAREMAP_H
#define DFMUXHARDWAREMAP_H

#include <G3Frame.h>
#include <G3Map.h>

/*
 * Class representing a single DfMux channel mapping. Module and channel numbers
 * are 0-indexed.
 */
class DfMuxChannelMapping : public G3FrameObject
{
public:
	int board_ip;
	int board_serial;
	int board_slot;
	int crate_serial;
	int module;
	int channel;

	std::string Description() const;
	std::string Summary() const;
	template <class A> void serialize(A &ar, unsigned v);
};

/* Map from bolometer ID to wiring info */
G3MAP_OF(std::string, DfMuxChannelMapping, DfMuxWiringMap);

G3_POINTERS(DfMuxChannelMapping);
G3_SERIALIZABLE(DfMuxChannelMapping, 2);

#endif

