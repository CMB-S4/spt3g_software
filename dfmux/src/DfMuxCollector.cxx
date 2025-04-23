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

#include <core/SetThreadName.h>
#include <core/G3Units.h>

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

	uint8_t num_modules; /* correct for V4 */
	uint8_t channels_per_module; /* sub-module block in V4 */
	uint8_t fir_stage;
	uint8_t module; /* linear; 0-7. don't like it much. */

	uint32_t seq; /* incrementing sequence number */

	int32_t s[256]; /* >= largest number of channels we expect */
	struct RawTimestamp ts;
} __attribute__((packed));

// Convert time stamp to a code in IRIG-B ticks (10 ns intervals)
static int64_t
RawTimestampToTimeCode(RawTimestamp stamp, double clock_scale)
{
	static __thread double last_code = -1;
	static __thread RawTimestamp last_stamp;
	static __thread double last_ss;
	struct tm tm;
	double ss;

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

	// fix the subsecond counter if the iceboard is using
	// an internal clock with a rate other than 100 MHz
	ss = le32toh(stamp.ss) * clock_scale;
	
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
		last_code = (last_code - last_ss) + ss;
	} else {
		tm.tm_mon = 0;       // Fake out timegm with the 274th of Jan.
		tm.tm_mday = tm.tm_yday; // since it ignores tm_yday otherwise
		last_code = 100000000LL * int64_t(timegm(&tm));
		last_code += ss;
	}

	last_stamp = stamp;
	last_ss = ss;
	
	return last_code;
}

DfMuxCollector::DfMuxCollector(G3EventBuilderPtr builder,
    std::vector<std::string> hosts) :
     builder_(builder), success_(false), stop_listening_(false)
{

	SetClockRate(100 * G3Units::MHz);
	success_ = (SetupSCTPSocket(hosts) != 0);
}

DfMuxCollector::DfMuxCollector(const char *listenaddr,
   G3EventBuilderPtr builder, std::vector<int32_t> board_list) :
    builder_(builder), success_(false), stop_listening_(false),
    board_list_(board_list)
{
	SetClockRate(100 * G3Units::MHz);
	success_ = (SetupUDPSocket(listenaddr) != 0);
}

DfMuxCollector::DfMuxCollector(const char *listenaddr,
    G3EventBuilderPtr builder, std::map<in_addr_t, int32_t> board_serial_map) :
     builder_(builder), success_(false), stop_listening_(false),
     board_serials_(board_serial_map)
{
	for (auto i : board_serials_)
		board_list_.push_back(i.second);

	SetClockRate(100 * G3Units::MHz);
	success_ = (SetupUDPSocket(listenaddr) != 0);
}

void DfMuxCollector::SetClockRate(double rate)
{

	clock_rate_ = rate;
	clock_scale_ = 100 * G3Units::MHz / clock_rate_;
}

int DfMuxCollector::SetupUDPSocket(const char *listenaddr)
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

int DfMuxCollector::SetupSCTPSocket(std::vector<std::string> hosts)
{

	fd_ = socket(AF_INET, SOCK_SEQPACKET, IPPROTO_SCTP);

	// Connect the socket to every board in our list. Errors are fatal.
	for (auto i : hosts) {
		struct addrinfo hints, *res;
		bzero(&hints, sizeof(hints));
		hints.ai_family = AF_INET;
		hints.ai_socktype = SOCK_SEQPACKET;
		hints.ai_protocol = 0;

		if (getaddrinfo(i.c_str(), "9876", &hints, &res))
			log_fatal("Could not resolve board \"%s\"", i.c_str());

		if (connect(fd_, res->ai_addr, res->ai_addrlen) != 0)
			log_fatal("Could not connect to board \"%s\" by SCTP "
			          "(%s). Maybe it has UDP-only firmware or is "
				  "not connected/powered?",
				  i.c_str(), strerror(errno));

		freeaddrinfo(res);
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
	setThreadName("dfmux listen");

	struct sockaddr_in addr;
	socklen_t addrlen = sizeof(addr);
	struct DfmuxPacket buf;
	ssize_t len;
	size_t target_size;
	int nchan;

	size_t base_size = sizeof(buf) - sizeof(buf.s);
	
	bzero(&addr, sizeof(addr));

	while (!collector->stop_listening_) {
		len = recvfrom(collector->fd_, &buf, sizeof(buf), 0,
		    (struct sockaddr *)&addr, &addrlen);
		nchan = (le16toh(buf.version) == 4) ? 128 : buf.channels_per_module;
		target_size = base_size + nchan*sizeof(buf.s[0])*2;
		if (len != (ssize_t)target_size) {
			log_error("Badly-sized packet with %d channels from %s "
			    "(%zd bytes should be %zd)",
			    nchan, inet_ntoa(addr.sin_addr),
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
	int mod_id;

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
	} else if (le16toh(packet->version) == 3 || le16toh(packet->version) == 4) {
		board_id = le16toh(packet->serial);
		if (board_list_.size() > 0 &&
		    std::find(board_list_.begin(), board_list_.end(), board_id)
		    == board_list_.end()) {
			struct in_addr i;
			i.s_addr = listenaddr_;
			log_debug("Received V%d data for board %d not "
			    "enumerated on listener for interface %s",
			    le16toh(packet->version),
			    board_id, inet_ntoa(i));
			return (-1);
		}
	} else {
		log_error("Unknown packet version %d from %s",
		    le16toh(packet->version), inet_ntoa(src));
		return (-1);
	}

	modseq = &sequence_[board_id];
	mod_id = (le16toh(packet->version) == 4) ?
	  (le32toh(packet->module) * 8 + le32toh(packet->channels_per_module)) :
	  le32toh(packet->module);
	seq = modseq->find(mod_id);
	if (seq != modseq->end()) {
		seq->second++;
		if ((uint32_t)seq->second != le32toh(packet->seq)) {
			log_warn("Out-of-order packet from %d/%d (%d "
			    "instead of %d)", board_id, mod_id,
			    le32toh(packet->seq), seq->second);
			seq->second = le32toh(packet->seq);
		}
	} else {
		// New board we haven't seen before
		(*modseq)[mod_id] = le32toh(packet->seq);
	}

	// Decode packet
	timecode = RawTimestampToTimeCode(packet->ts, clock_scale_);

	// All times reported by the readout are exactly one second behind,
	// likely due to a misparsing of which IRIG code the time marker refers
	// to (before or after). Shift them all 1 second forward.
	timecode += 100000000LL;

	int nchan = (le16toh(packet->version) == 4) ? 128 : le32toh(packet->channels_per_module);
	DfMuxSamplePtr sample(new DfMuxSample(timecode, nchan*2));

	// NB: Bottom 8 bits are zero. Divide by 256 rather than >> 8 to
	// guarantee sign preservation.
	for (int i = 0; i < sample->NSamples(); i++)
		sample->Samples()[i] = int32_t(le32toh(packet->s[i])) / 256;

	DfMuxSamplePacketPtr outpacket(new DfMuxSamplePacket);
	outpacket->board = board_id;
	outpacket->module = le32toh(packet->module);
	outpacket->block = (le16toh(packet->version) == 4) ?
	  le32toh(packet->channels_per_module) : 0;
	outpacket->sample = sample;

	// Work around bug in IceBoard firmware. You always get 8 modules'
	// worth of packets, no matter what packet->num_modules says.
	// (firmware V3 and below)
	outpacket->nmodules = (le16toh(packet->version) == 4) ? le32toh(packet->num_modules) : 8;
	outpacket->nblocks = (le16toh(packet->version) == 4) ? 8 : 1;
	outpacket->nchannels = nchan;

	builder_->AsyncDatum(timecode, outpacket);

	return 0;
}

static DfMuxCollectorPtr
make_dfmux_collector_v2_from_dict(const char *listenaddr,
    G3EventBuilderPtr builder, py::dict board_serial_map)
{
	std::map<in_addr_t, int32_t> serial_map;
	
	py::list items = board_serial_map.items();
	for (ssize_t i = 0; i < py::len(items); i++) {
		int32_t serial = py::extract<int>(items[i][1])();
		in_addr_t ip;
		bool found = false;

		auto int_item = py::extract<int>(items[i][0]);
		if (int_item.check()) {
			ip = int_item();
			found = true;
		} else {
			auto str_item = py::extract<std::string>(items[i][0]);
			if (str_item.check()) {
				std::string host = str_item();
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

				found = true;
			}
		}
		if (!found) {
			log_fatal("Map keys must be integer or string "
			    "representations of the IP address or hostname");
		}

		serial_map[ip] = serial;
	}

	return DfMuxCollectorPtr(new DfMuxCollector(listenaddr, builder,
	    serial_map));
}

PYBINDINGS("dfmux", scope)
{
	register_class<DfMuxCollector>(scope, "DfMuxCollector",
	    "Listener object that collects IceBoard packets from a single network "
	    "interface and decodes and forwards them to a DfMuxBuilder object for "
	    "insertion into the data stream. Takes the builder object to which the "
	    "data should be sent. Uses either multicast UDP or SCTP depending on "
	    "which constructor is used.")
	    .def(py::init<G3EventBuilderPtr, std::vector<std::string> >(),
	        py::arg("builder"), py::arg("hostnames"),
	        "Create a DfMuxCollector listening for SCTP packets from the "
	        "listed hosts (e.g. [\"iceboard0062.local\", ...]) and forwards "
	        "it to DfMuxBuilder \"builder\".")
	    .def(py::init<const char *, G3EventBuilderPtr, std::vector<int32_t> >(),
	        py::arg("interface"), py::arg("builder"),
	        py::arg("boardlist")=std::vector<int32_t>(),
	        "Create a DfMuxCollector listening on \"interface\" for multicasted "
	        "UDP packets and forwards it to DfMuxBuilder \"builder\". Filters to "
	        "only the boards specified in \"boardlist\" (by default empty, "
	        "implying all boards).")
	    .def("__init__", py::make_constructor(make_dfmux_collector_v2_from_dict,
	        py::default_call_policies(),
	        (py::arg("interface"), py::arg("builder"), py::arg("board_serial_map"))),
	        "Crate a DfMuxCollector that can parse V2 (64x) data. Pass a mapping "
	        "from board IP address (strings or integers) to serial numbers as "
	        "the last argument")
	    .def("Start", &DfMuxCollector::Start)
	    .def("Stop", &DfMuxCollector::Stop)
	    .def_property("clock_rate", &DfMuxCollector::GetClockRate,
	        &DfMuxCollector::SetClockRate, "Set the clock rate for the "
	        "iceboard subseconds counter, e.g. for hidfmux.  Values should be "
	        "in G3Units of frequency.  Defaults to 100*core.G3Units.MHz.")
	;
}

