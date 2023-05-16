#ifndef DFMUX_DATA_COLLECTOR_H
#define DFMUX_DATA_COLLECTOR_H

#include <map>
#include <arpa/inet.h>
#include <G3EventBuilder.h>

#include <dfmux/DfMuxSample.h>

/*
 * DfMuxSamplePacket: Internal class used to communicate data from
 * DfMuxCollector to DfMuxBuilder. Does not appear in frames.
 */

class DfMuxSamplePacket : public G3FrameObject {
public:
	int32_t board;		/* Board serial number */
	int32_t module;		/* Module number (0-7) */
	int32_t nmodules;	/* Total number of modules on board (8) */
	DfMuxSamplePtr sample;	/* Pointer to the DfMuxSample */
};

G3_POINTERS(DfMuxSamplePacket);


/*
 * DfMuxCollector: Listens to multicast data from ICE board streamers
 * and forwards to a DfMuxBuilder. 
 */

class DfMuxCollector {
public:
	// Listen for SCTP data arriving from the hostnames/IPs given
	// in the boards_to_listen_for argument.
	DfMuxCollector(G3EventBuilderPtr builder,
	    std::vector<std::string> boards_to_listen_for);

	// Listen for multicast UDP data arriving on listenaddr from the
	// boards listed. If no boards are listed, will forward data
	// from all boards transmitting to listenaddr.
	DfMuxCollector(const char *listenaddr, G3EventBuilderPtr builder,
	    std::vector<int32_t> boards_to_listen_for);

	// This constructor takes an IP address -> board serial mapping,
	// which is required to process V2 (64x) streamer packets.
	DfMuxCollector(const char *listenaddr, G3EventBuilderPtr builder,
	    std::map<in_addr_t, int32_t> board_serial_map);
	~DfMuxCollector();

	int Start();	// Start listening thread
	int Stop();	// Stop listening thread

	// Configure iceboard clock rate
	void SetClockRate(uint64_t rate) { clock_rate_ = rate; };
	uint64_t GetClockRate() { return clock_rate_; };

private:
	int SetupUDPSocket(const char *listenaddr);
	int SetupSCTPSocket(std::vector<std::string> hosts);
	static void Listen(DfMuxCollector *collector);
	std::thread listen_thread_;

	int BookPacket(struct DfmuxPacket *packet, struct in_addr src);

	std::map<int32_t, std::map<int32_t, int32_t> > sequence_;

	G3EventBuilderPtr builder_;
	bool success_;
	volatile bool stop_listening_;
	std::map<in_addr_t, int32_t> board_serials_;
	std::vector<int32_t> board_list_;
	int fd_;
	in_addr_t listenaddr_;
	uint64_t clock_rate_;

	SET_LOGGER("DfMuxCollector");
};

G3_POINTERS(DfMuxCollector);

#endif // DFMUX_DATA_COLLECTOR_H

