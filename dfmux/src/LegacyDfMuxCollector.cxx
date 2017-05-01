#include <pybindings.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <boost/python.hpp>

#ifdef __FreeBSD__
#include <sys/endian.h>
#endif

#ifdef __APPLE__
	
#include <libkern/OSByteOrder.h>

#define htobe16(x) OSSwapHostToBigInt16(x)
#define htole16(x) OSSwapHostToLittleInt16(x)
#define be16toh(x) OSSwapBigToHostInt16(x)
#define le16toh(x) OSSwapLittleToHostInt16(x)

#define htobe32(x) OSSwapHostToBigInt32(x)
#define htole32(x) OSSwapHostToLittleInt32(x)
#define be32toh(x) OSSwapBigToHostInt32(x)
#define le32toh(x) OSSwapLittleToHostInt32(x)

#define htobe64(x) OSSwapHostToBigInt64(x)
#define htole64(x) OSSwapHostToLittleInt64(x)
#define be64toh(x) OSSwapBigToHostInt64(x)
#define le64toh(x) OSSwapLittleToHostInt64(x)
#endif

#include <dfmux/DfMuxCollector.h>

union RawTimestamp {
	/* SIGNED types are used here to permit negative numbers during
	 * renormalization and comparison. THE ORDER IS IMPORTANT, since
	 * this matches the VHDL's field packing. */
	struct {
		int32_t s;
		int32_t ss;
	} irig_test;
	struct {
		int32_t y,d,h,m,s;
		int32_t ss;
	} irig;
	struct {
		uint32_t source;
		int32_t maj;
		int32_t ticks;
		uint32_t recent;
	} timesync;
};

#define MUX_MAGIC		0x666d7578
#define MUX_VERSION		1
#define STREAMER_IRIG		1
#define STREAMER_IRIG_TEST	2
#define STREAMER_TIMESYNC	3

struct DfmuxPacket {
	uint32_t magic;
	uint32_t version;
	uint32_t reserved[2];

	int32_t s[128];

	uint32_t ts_port;
	union RawTimestamp ts;
} __attribute__((packed));

// Convert time stamp to a code in IRIG-B ticks (10 ns intervals)
static int64_t
RawTimestampToTimeCode(RawTimestamp stamp, uint32_t port)
{
	static __thread int64_t last_code = -1;
	static __thread RawTimestamp last_stamp;
	struct tm tm;

	if (port == STREAMER_IRIG_TEST) {
		int64_t code;
		code = 100000000LL * int64_t(be32toh(stamp.irig_test.s));
		code += (uint64_t)be32toh(stamp.irig_test.ss);

		return code;
	}
	
	g3_assert(port == STREAMER_IRIG);

	// Some IRIG generators don't fill in year, which is annoying. This
	// algorithm works unless time goes backwards, the computer clock's year
	// is wrong, or data acquisition begins within a ms of the new year.
	if (be32toh(stamp.irig.y) == 0) {
		if (last_code == -1) {
			// Assume the host computer at least has the year right
			time_t curtime;
			curtime = time(NULL);
			gmtime_r(&curtime, &tm);

			stamp.irig.y = htobe32(tm.tm_year % 100);
		} else {
			// Track from the last stamp. If New Year's seems to
			// have happened, increment the year.

			if (be32toh(stamp.irig.d) == 1 &&
			    be32toh(last_stamp.irig.d) >= 365)
				stamp.irig.y = htobe32(be32toh(
				    last_stamp.irig.y) + 1);
			else
				stamp.irig.y = last_stamp.irig.y;
		}
	}
	
	tm.tm_year = be32toh(stamp.irig.y) + 100 /* tm_year starts in 1900 */;
	tm.tm_yday = be32toh(stamp.irig.d);
	tm.tm_hour = be32toh(stamp.irig.h);
	tm.tm_min = be32toh(stamp.irig.m);
	tm.tm_sec = be32toh(stamp.irig.s);

	if (last_code != -1 && stamp.irig.y == last_stamp.irig.y &&
	    stamp.irig.d == last_stamp.irig.d &&
	    stamp.irig.h == last_stamp.irig.h &&
	    stamp.irig.m == last_stamp.irig.m &&
	    stamp.irig.s == last_stamp.irig.s) {
		// If all fields but sub-second agree, just apply the change
		// in the subsecond field as an offset
		last_code = (last_code - be32toh(last_stamp.irig.ss)) +
		    be32toh(stamp.irig.ss);
	} else {
		tm.tm_mon = 0;       // Fake out timegm with the 274th of Jan.
		tm.tm_mday = tm.tm_yday; // since it ignores tm_yday otherwise
		last_code = 100000000LL * int64_t(timegm(&tm));
		last_code += (uint64_t)be32toh(stamp.irig.ss);
	}

	last_stamp = stamp;
	
	return last_code;
}

class LegacyDfMuxCollector {
public:
        LegacyDfMuxCollector(int port, G3EventBuilderPtr builder,
	    const char *mcast_listen_addr = "", const char *mcast_group = "");
        ~LegacyDfMuxCollector();

        int Start();
        int Stop();

private:
        static void Listen(LegacyDfMuxCollector *collector);
        std::thread listen_thread_;

        int BookPacket(struct DfmuxPacket *packet, struct in_addr src);

        G3EventBuilderPtr builder_;
        bool success_;
        volatile bool stop_listening_;
        int fd_;

        SET_LOGGER("LegacyDfMuxCollector");
};

LegacyDfMuxCollector::LegacyDfMuxCollector(int port, G3EventBuilderPtr builder,
    const char *mcast_listen_addr, const char *mcast_group)
    : builder_(builder), success_(false), stop_listening_(false)
{
	struct sockaddr_in addr;
	struct ip_mreq mcast;
	int yes;

	fd_ = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);

	// Allow multiple listeners
	yes = 1;
#ifdef __linux__
	if (setsockopt(fd_, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes)) < 0)
		perror("Error setting SO_REUSEADDR");
#else
	if (setsockopt(fd_, SOL_SOCKET, SO_REUSEPORT, &yes, sizeof(yes)) < 0)
		perror("Error setting SO_REUSEPORT");
#endif

	addr.sin_family = AF_INET;
	addr.sin_addr.s_addr = htonl(INADDR_ANY);
	addr.sin_port = htons(port);
	if (bind(fd_, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
		perror(NULL);
		return;
	}

	if (mcast_group != NULL && mcast_group[0] != '\0') {
		mcast.imr_multiaddr.s_addr = inet_addr(mcast_group);
		mcast.imr_interface.s_addr = inet_addr(mcast_listen_addr);

		if (setsockopt(fd_, IPPROTO_IP, IP_ADD_MEMBERSHIP, &mcast,
		    sizeof(mcast)) < 0) {
			perror(NULL);
			return;
		}
	}

	int rcvbuf = 80000 * sizeof(struct DfmuxPacket);
	if (setsockopt(fd_, SOL_SOCKET, SO_RCVBUF, &rcvbuf, sizeof(rcvbuf)) < 0)
		perror("Error setting receive queue length");

	success_ = true;
}

LegacyDfMuxCollector::~LegacyDfMuxCollector()
{
	Stop();
	close(fd_);
}

int LegacyDfMuxCollector::Start()
{
	stop_listening_ = false;
	listen_thread_ = std::thread(Listen, this);

	return (0);
}

int LegacyDfMuxCollector::Stop()
{
	stop_listening_ = true;
	listen_thread_.join();

	return (0);
}

void LegacyDfMuxCollector::Listen(LegacyDfMuxCollector *collector)
{
	struct sockaddr_in addr;
	socklen_t addrlen = sizeof(addr);
	struct DfmuxPacket buf;
	ssize_t len;
	
	bzero(&addr, sizeof(addr));

	while (!collector->stop_listening_) {
		len = recvfrom(collector->fd_, &buf, sizeof(buf), 0,
		    (struct sockaddr *)&addr, &addrlen);
		if (len != sizeof(buf)) {
			log_error("Badly-sized packet from %s "
			    "(%zd bytes should be %zd)",
			    inet_ntoa(addr.sin_addr), len, sizeof(buf));
			continue;
		}

		collector->BookPacket(&buf, addr.sin_addr);
	}
}

int LegacyDfMuxCollector::BookPacket(struct DfmuxPacket *packet,
    struct in_addr src)
{
	int64_t timecode;

	if (be32toh(packet->magic) != MUX_MAGIC) {
		log_error("Corrupted packet from %s begins with %#x "
		    "instead of %#x", inet_ntoa(src), be32toh(packet->magic),
		    MUX_MAGIC);
		return (-1);
	}

	// Decode packet
	timecode = RawTimestampToTimeCode(packet->ts, be32toh(packet->ts_port));

	// Split into four apparent packets, one for each module.
	// This matches the convention for the ICE system

	for (int m = 0; m < 4; m++) {
		DfMuxSamplePtr sample(new DfMuxSample(timecode, 32));

		// NB: Bottom 8 bits are zero. Divide by 256 rather than >> 8
		// to guarantee sign preservation.
		for (int i = 0; i < sample->NSamples(); i++)
			sample->Samples()[i] =
			    int32_t(be32toh(packet->s[i + m*32])) / 256;

		DfMuxSamplePacketPtr outpacket(new DfMuxSamplePacket);

		// Use last octet of IP address as "serial number"
		outpacket->board = ntohl(src.s_addr) & 0xff;
	
		outpacket->sample = sample;
		outpacket->module = m;
		outpacket->nmodules = 4;

		builder_->AsyncDatum(timecode, outpacket);
	}

	return 0;
}

PYBINDINGS("dfmux")
{
	namespace bp = boost::python;

	bp::class_<LegacyDfMuxCollector, boost::shared_ptr<LegacyDfMuxCollector>, boost::noncopyable>("LegacyDfMuxCollector", "Listener object that collects legacy network packets and decodes and forwards them to a DfMuxBuilder object for insertion into the data stream. Takes a builder object to which the data should be sent. If boards running in multicast code, also takes the address of the interface on which to listen and the multicast group address. Only works with DAN streamers.", bp::init<int, G3EventBuilderPtr, const char *, const char *>((bp::arg("port"), bp::arg("builder"), bp::arg("mcastlistenaddr")="", bp::arg("mcastgroupaddr")="")))
	    .def("Start", &LegacyDfMuxCollector::Start)
	    .def("Stop", &LegacyDfMuxCollector::Stop)
	;
}

