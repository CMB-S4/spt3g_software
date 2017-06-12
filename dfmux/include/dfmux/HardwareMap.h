#ifndef DFMUXHARDWAREMAP_H
#define DFMUXHARDWAREMAP_H

#include <G3Frame.h>
#include <G3Map.h>

#include <stdint.h>

/*
 * Class representing a single DfMux channel mapping. Module and channel numbers
 * are 0-indexed.
 */
class DfMuxChannelMapping : public G3FrameObject
{
public:
	int32_t board_ip;
	int32_t board_serial;
	int32_t board_slot;
	int32_t crate_serial;
	int32_t module;
	int32_t channel;

	std::string Description() const;
	std::string Summary() const;
	template <class A> void serialize(A &ar, unsigned v);
};

/* Map from bolometer ID to wiring info */
G3MAP_OF(std::string, DfMuxChannelMapping, DfMuxWiringMap);

G3_POINTERS(DfMuxChannelMapping);
G3_SERIALIZABLE(DfMuxChannelMapping, 2);

#endif

