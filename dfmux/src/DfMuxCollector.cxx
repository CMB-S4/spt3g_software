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
#include <netdb.h>
#include <boost/python.hpp>

#ifdef __FreeBSD__
#include <sys/endian.h>
#include <pthread_np.h>
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


struct RawTimestamp {
	/* SIGNED types are used here to permit negative numbers during
	 * renormalization and comparison. THE ORDER IS IMPORTANT, since
	 * this matches the VHDL's field packing. */
	int32_t y,d,h,m,s;
	int32_t ss;
	int32_t c, sbs;
};

#define FAST_MAGIC	0x666d7578
struct DfmuxPacket {
	uint32_t magic;
	uint16_t version;

	uint16_t serial; /* zero for V2, filled for V3 */

	uint8_t num_modules;
	uint8_t channels_per_module;
	uint8_t fir_stage;
	uint8_t module; /* linear; 0-7. don't like it much. */

	uint32_t seq; /* incrementing sequence number */

	int32_t s[256]; /* >= largest number of channels we expect */
	struct RawTimestamp ts;
} __attribute__((packed));

// Convert time stamp to a code in IRIG-B ticks (10 ns intervals)
static int64_t
RawTimestampToTimeCode(RawTimestamp stamp)
{
	static __thread int64_t last_code = -1;
	static __thread RawTimestamp last_stamp;
	struct tm tm;

	// Some IRIG generators don't fill in year, which is annoying. This
	// algorithm works unless time goes backwards, the computer clock's year
	// is wrong, or data acquisition begins within a ms of the new year.
	if (le32toh(stamp.y) == 0) {
		if (last_code == -1) {
			// Assume the host computer at least has the year right
			time_t curtime;
			curtime = time(NULL);
			gmtime_r(&curtime, &tm);

			stamp.y = htole32(tm.tm_year % 100);
		} else {
			// Track from the last stamp. If New Year's seems to
			// have happened, increment the year.

			if (le32toh(stamp.d) == 1 &&
			    le32toh(last_stamp.d) >= 365)
				stamp.y = htole32(le32toh(last_stamp.y) + 1);
			else
				stamp.y = last_stamp.y;
		}
	}
	
	tm.tm_year = le32toh(stamp.y) + 100 /* tm_year starts in 1900 */;
	tm.tm_yday = le32toh(stamp.d);
	tm.tm_hour = le32toh(stamp.h);
	tm.tm_min = le32toh(stamp.m);
	tm.tm_sec = le32toh(stamp.s);

	if (last_code != -1 && stamp.y == last_stamp.y &&
	    stamp.d == last_stamp.d && stamp.h == last_stamp.h &&
	    stamp.m == last_stamp.m && stamp.s == last_stamp.s) {
		// If all fields but sub-second agree, just apply the change
		// in the subsecond field as an offset
		last_code = (last_code - le32toh(last_stamp.ss)) +
		    le32toh(stamp.ss);
	} else {
		tm.tm_mon = 0;       // Fake out timegm with the 274th of Jan.
		tm.tm_mday = tm.tm_yday; // since it ignores tm_yday otherwise
		last_code = 100000000LL * int64_t(timegm(&tm));
		last_code += (uint64_t)le32toh(stamp.ss);
	}

	last_stamp = stamp;
	
	return last_code;
}

DfMuxCollector::DfMuxCollector(const char *listenaddr,
   G3EventBuilderPtr builder, std::vector<int32_t> board_list) :
    builder_(builder), success_(false), stop_listening_(false),
    board_list_(board_list)
{
	success_ = (SetupSocket(listenaddr) != 0);
}

DfMuxCollector::DfMuxCollector(const char *listenaddr,
    G3EventBuilderPtr builder, std::map<in_addr_t, int32_t> board_serial_map) :
     builder_(builder), success_(false), stop_listening_(false),
     board_serials_(board_serial_map)
{
	for (auto i : board_serials_)
		board_list_.push_back(i.second);

	success_ = (SetupSocket(listenaddr) != 0);
}

int DfMuxCollector::SetupSocket(const char *listenaddr)
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

	// Listen on 239.192.0.2, port 9876
	addr.sin_family = AF_INET;
	addr.sin_addr.s_addr = htonl(INADDR_ANY);
	addr.sin_port = htons(9876);
	if (bind(fd_, (struct sockaddr *)&addr, sizeof(addr)) < 0) {
		perror(NULL);
		return errno;
	}

	mcast.imr_multiaddr.s_addr = htonl(0xefc00002 /* 239.192.0.2 */);
	listenaddr_ = inet_addr(listenaddr);
	mcast.imr_interface.s_addr = listenaddr_;

	if (setsockopt(fd_, IPPROTO_IP, IP_ADD_MEMBERSHIP, &mcast,
	    sizeof(mcast)) < 0) {
		perror(NULL);
		return errno;
	}

	int rcvbuf = 80000 * sizeof(struct DfmuxPacket);
	if (setsockopt(fd_, SOL_SOCKET, SO_RCVBUF, &rcvbuf, sizeof(rcvbuf)) < 0)
		perror("Error setting receive queue length");

	return 0;
}

DfMuxCollector::~DfMuxCollector()
{
	Stop();
	close(fd_);
}

int DfMuxCollector::Start()
{
	stop_listening_ = false;
	listen_thread_ = std::thread(Listen, this);

#ifdef __linux__
	pthread_setname_np(listen_thread_.native_handle(), "dfmux listen");
#endif
#ifdef __FreeBSD__
	pthread_set_name_np(listen_thread_.native_handle(), "dfmux listen");
#endif

	return (0);
}

int DfMuxCollector::Stop()
{
	stop_listening_ = true;
	listen_thread_.join();

	return (0);
}

void DfMuxCollector::Listen(DfMuxCollector *collector)
{
	struct sockaddr_in addr;
	socklen_t addrlen = sizeof(addr);
	struct DfmuxPacket buf;
	ssize_t len;
	size_t target_size;

	size_t base_size = sizeof(buf) - sizeof(buf.s);
	
	bzero(&addr, sizeof(addr));

	while (!collector->stop_listening_) {
		len = recvfrom(collector->fd_, &buf, sizeof(buf), 0,
		    (struct sockaddr *)&addr, &addrlen);
		target_size = base_size +
		    buf.channels_per_module*sizeof(buf.s[0])*2;
		if (len != target_size) {
			log_error("Badly-sized packet with %d channels from %s "
			    "(%zd bytes should be %zd)",
			    buf.channels_per_module, inet_ntoa(addr.sin_addr),
			    len, target_size);
			continue;
		}
		
		// Move timestamp from packed location to correct field
		memmove(&buf.ts, (char *)&buf.ts - (sizeof(buf) - target_size),
		    sizeof(buf.ts));

		collector->BookPacket(&buf, addr.sin_addr);
	}
}

int DfMuxCollector::BookPacket(struct DfmuxPacket *packet, struct in_addr src)
{
	std::map<int32_t, int32_t> *modseq;
	std::map<int32_t, int32_t>::iterator seq;
	int64_t timecode;
	int board_id;

	if (le32toh(packet->magic) != FAST_MAGIC) {
		log_error("Corrupted packet from %s begins with %#x "
		    "instead of %#x", inet_ntoa(src), le32toh(packet->magic),
		    FAST_MAGIC);
		return (-1);
	}

	if (le16toh(packet->version) == 2) {
		auto id_it = board_serials_.find(src.s_addr);
		if (id_it == board_serials_.end()) {
			log_debug("V2 (64x) data from unknown board %s.",
			    inet_ntoa(src));
			return (-1);
		}
		board_id = id_it->second;
	} else if (le16toh(packet->version) == 3) {
		board_id = le16toh(packet->serial);
		if (board_list_.size() > 0 &&
		    std::find(board_list_.begin(), board_list_.end(), board_id)
		    == board_list_.end()) {
			struct in_addr i;
			i.s_addr = listenaddr_;
			log_debug("Received V3 data for board %d not "
			    "enumerated on listener for interface %s",
			    board_id, inet_ntoa(i));
			return (-1);
		}
	} else {
		log_error("Unknown packet version %d from %s",
		    le16toh(packet->version), inet_ntoa(src));
		return (-1);
	}

	modseq = &sequence_[board_id];
	seq = modseq->find(le32toh(packet->module));
	if (seq != modseq->end()) {
		seq->second++;
		if ((uint32_t)seq->second != le32toh(packet->seq)) {
			log_warn("Out-of-order packet from %d/%d (%d "
			    "instead of %d)", board_id, le32toh(packet->module),
			    le32toh(packet->seq), seq->second);
			seq->second = le32toh(packet->seq);
		}
	} else {
		// New board we haven't seen before
		(*modseq)[le32toh(packet->module)] = le32toh(packet->seq); 
	}

	// Decode packet
	timecode = RawTimestampToTimeCode(packet->ts);

	// All times reported by the readout are exactly one second behind,
	// likely due to a misparsing of which IRIG code the time marker refers
	// to (before or after). Shift them all 1 second forward.
	timecode += 100000000LL;

	DfMuxSamplePtr sample(new DfMuxSample(timecode,
	    le32toh(packet->channels_per_module)*2));

	// NB: Bottom 8 bits are zero. Divide by 256 rather than >> 8 to
	// guarantee sign preservation.
	for (int i = 0; i < sample->NSamples(); i++)
		sample->Samples()[i] = int32_t(le32toh(packet->s[i])) / 256;

	DfMuxSamplePacketPtr outpacket(new DfMuxSamplePacket);
	outpacket->board = board_id;
	outpacket->module = le32toh(packet->module);
	outpacket->sample = sample;

	// Work around bug in IceBoard firmware. You always get 8 modules'
	// worth of packets, no matter what packet->num_modules says.
	outpacket->nmodules = 8;

	builder_->AsyncDatum(timecode, outpacket);

	return 0;
}

static DfMuxCollectorPtr
make_dfmux_collector_v2_from_dict(const char *listenaddr,
    G3EventBuilderPtr builder, boost::python::dict board_serial_map)
{
	namespace bp = boost::python;
	std::map<in_addr_t, int32_t> serial_map;
	
	bp::list items = board_serial_map.items();
	for (size_t i = 0; i < bp::len(items); i++) {
		int32_t serial = bp::extract<int>(items[i][1])();
		in_addr_t ip;

		if (bp::extract<int>(items[i][0]).check()) {
			ip = bp::extract<int>(items[i][0])();
		} else if (bp::extract<std::string>(items[i][0]).check()) {
			std::string host = bp::extract<std::string>(
			    items[i][0])();
			struct addrinfo hints, *info;
			int err;

			bzero(&hints, sizeof(hints));
			hints.ai_family = PF_INET;
			
			err = getaddrinfo(host.c_str(), NULL, &hints, &info);
			if (err != 0)
				log_fatal("Could not find host %s (%s)",
				    host.c_str(), gai_strerror(err));

			g3_assert(info->ai_family == PF_INET);
			ip = ((struct sockaddr_in *)(info->ai_addr))->
			    sin_addr.s_addr;
		} else {
			log_fatal("Map keys must be integer or string "
			    "representations of the IP address or hostname");
		}

		serial_map[ip] = serial;
	}

	return DfMuxCollectorPtr(new DfMuxCollector(listenaddr, builder,
	    serial_map));
}

PYBINDINGS("dfmux")
{
	namespace bp = boost::python;

	bp::class_<DfMuxCollector, DfMuxCollectorPtr, boost::noncopyable>("DfMuxCollector", "Listener object that collects IceBoard packets from a single network interface and decodes and forwards them to a DfMuxBuilder object for insertion into the data stream. Takes the IP address of the interface on the host computer on which to listen as well as the builder object to which the data should be sent.", bp::init<const char *, G3EventBuilderPtr, std::vector<int32_t> >((bp::arg("interface"), bp::arg("builder"), bp::arg("boardlist")=std::vector<int32_t>()), "Create a DfMuxCollector listening on \"interface\" for multicasted packets and forwards it to DfMuxBuilder \"builder\". Filters to only the boards specified in \"boardlist\" (by default empty, implying all boards)."))
	    .def("__init__", bp::make_constructor(make_dfmux_collector_v2_from_dict, bp::default_call_policies(), (bp::arg("interface"), bp::arg("builder"), bp::arg("board_serial_map"))), "Crate a DfMuxCollector that can parse V2 (64x) data. Pass a mapping from board IP address (strings or integers) to serial numbers as the last argument")
	    .def("Start", &DfMuxCollector::Start)
	    .def("Stop", &DfMuxCollector::Stop)
	;
}

