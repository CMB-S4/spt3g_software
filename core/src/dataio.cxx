#include <G3Logging.h>
#include <dataio.h>
#include "streams.h"

#include <filesystem>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdlib.h>


static int
connect_remote(const std::string &path, float timeout)
{
	// TCP Socket. Two syntaxes:
	// - tcp://host:port -> connect to "host" on "port" and read
	//   until EOF
	// - tcp://*:port -> listen on "port" for the first connection
	//   and read until EOF

	std::string host = path.substr(path.find("://") + 3);
	if (host.find(":") == host.npos)
		log_fatal("Could not open URL %s: unspecified port", path.c_str());
	std::string port = host.substr(host.find(":") + 1);
	host = host.substr(0, host.find(":"));

	log_debug("Opening connection to %s, port %s", host.c_str(), port.c_str());

	int fd = -1;

	if (strcmp(host.c_str(), "*") == 0) {
		// Listen for incoming connections
		struct sockaddr_in6 sin;
		int no = 0, yes = 1;
		int lfd;

		bzero(&sin, sizeof(sin));
		sin.sin6_family = AF_INET6;
		sin.sin6_port = htons(strtol(port.c_str(), NULL, 10));
	#ifdef SIN6_LEN
		sin.sin6_len = sizeof(sin);
	#endif

		lfd = socket(PF_INET6, SOCK_STREAM, 0);
		if (lfd <= 0)
			log_fatal("Could not listen on %s (%s)",
			    path.c_str(), strerror(errno));
		setsockopt(lfd, IPPROTO_IPV6, IPV6_V6ONLY, &no, sizeof(no));
		setsockopt(lfd, SOL_SOCKET, SO_REUSEADDR, &yes, sizeof(yes));

		if (bind(lfd, (struct sockaddr *)&sin, sizeof(sin)) < 0)
			log_fatal("Could not bind on port %s (%s)",
			    port.c_str(), strerror(errno));
		if (listen(lfd, 1) < 0)
			log_fatal("Could not listen on port %s (%s)",
			    port.c_str(), strerror(errno));

		log_debug("Waiting for connection on port %s", port.c_str());
		fd = accept(lfd, NULL, NULL);
		log_debug("Accepted connection on port %s", port.c_str());
		close(lfd);

		return fd;
	}

	// Connect to a listening host elsewhere

	struct addrinfo hints, *info, *r;
	int err;

	bzero(&hints, sizeof(hints));
	hints.ai_family = AF_UNSPEC;
	hints.ai_socktype = SOCK_STREAM;

	err = getaddrinfo(host.c_str(), port.c_str(), &hints, &info);
	if (err != 0)
		log_fatal("Could not find host %s (%s)",
		    host.c_str(), gai_strerror(err));

	// Loop through possible addresses until we find one
	// that works.
	fd = -1;
	for (r = info; r != NULL; r = r->ai_next) {
		fd = socket(r->ai_family, r->ai_socktype, r->ai_protocol);
		if (fd == -1)
			continue;

		if (connect(fd, r->ai_addr, r->ai_addrlen) == -1) {
			close(fd);
			fd = -1;
			continue;
		}

		break;
	}

	if (fd == -1)
		log_fatal("Could not connect to %s (%s)",
		    path.c_str(), strerror(errno));

        if (timeout >= 0) {
                struct timeval tv;
                tv.tv_sec = (int)timeout;
                tv.tv_usec = (int)(1e6 * (timeout - tv.tv_sec));
                if (setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO,
                               (char *)&tv, sizeof(tv)) < 0)
                        log_fatal("Failed to set timeout on socket; errno=%i",
                                  errno);
        }

	if (info != NULL)
		freeaddrinfo(info);

	return fd;
}


enum Codec {
  NONE = 0,
  GZ = 1,
  BZIP2 = 2,
  LZMA = 3,
  REMOTE = 4,
};

static bool
has_ext(const std::string &path, const std::string &ext)
{
	size_t n = ext.size();
	if (path.size() <= n)
		return false;
	return !path.compare(path.size() - n, n, ext);
}

static Codec
get_codec(const std::string &path, const std::string &ext=".g3")
{
	if (has_ext(path, ext + ".gz")) {
#ifdef ZLIB_FOUND
		return GZ;
#else
		log_fatal("Library not compiled with gzip support");
#endif
	}
	if (has_ext(path, ext + ".bz2")) {
#ifdef BZIP2_FOUND
		return BZIP2;
#else
		log_fatal("Library not compiled with bzip2 support");
#endif
	}
	if (has_ext(path, ext + ".xz")) {
#ifdef LZMA_FOUND
		return LZMA;
#else
		log_fatal("Library not compiled with LZMA support");
#endif
	}

	if (!ext.size() || has_ext(path, ext))
		return NONE;

	log_fatal("Filename %s does not have extension %s",
	    path.c_str(), ext.c_str());
}

static Codec
check_input_path(const std::string &path, const std::string &ext)
{
	if (path.find("tcp://") == 0)
		return REMOTE;

	std::filesystem::path fpath(path);
	if (!std::filesystem::exists(fpath) ||
	    !std::filesystem::is_regular_file(fpath))
		log_fatal("Could not find file %s", path.c_str());

	return get_codec(path, ext);
}

static Codec
check_output_path(const std::string &path, const std::string &ext)
{
	std::filesystem::path fpath(path);

	if (fpath.empty())
		log_fatal("Empty file path");

	if (fpath.has_parent_path()) {
		auto ppath = fpath.parent_path();
		if (!std::filesystem::exists(ppath))
			log_fatal("Parent path does not exist: %s",
			    ppath.string().c_str());
		if (!std::filesystem::is_directory(ppath))
			log_fatal("Parent path is not a directory: %s",
			    ppath.string().c_str());
	}

	return get_codec(path, ext);
}


static void
reset_stream(std::ios &stream)
{
	std::streambuf* sbuf = stream.rdbuf();
	if (sbuf) {
		sbuf->pubsync();
		delete sbuf;
	}
	stream.rdbuf(nullptr);
	stream.pword(0) = nullptr;
}

static void
stream_cb(std::ios::event ev, std::ios_base& stream, int index)
{
	std::streambuf* buf = nullptr;

	switch (ev) {
	case std::ios::event::erase_event:
		buf = static_cast<std::streambuf*>(stream.pword(0));
		if (buf) {
			buf->pubsync();
			delete buf;
			stream.pword(0) = nullptr;
		}
		break;
	default:
		break;
	}
}

void
g3_istream_from_path(std::istream &stream, const std::string &path, float timeout,
    size_t buffersize, const std::string &ext)
{
	reset_stream(stream);

	int fd = -1;

	// Figure out what kind of ultimate data source this is
	switch(check_input_path(path, ext)) {
	case REMOTE:
		fd = connect_remote(path, timeout);
		stream.rdbuf(new RemoteInputStreamBuffer(fd, buffersize));
		break;
#ifdef ZLIB_FOUND
	case GZ:
		stream.rdbuf(new GZipDecoder(path, buffersize));
		break;
#endif
#ifdef BZIP2_FOUND
	case BZIP2:
		stream.rdbuf(new BZip2Decoder(path, buffersize));
		break;
#endif
#ifdef LZMA_FOUND
	case LZMA:
		stream.rdbuf(new LZMADecoder(path, buffersize));
		break;
#endif
	default:
		// Read buffer
		stream.rdbuf(new InputFileStreamCounter(path, buffersize));
		break;
	}

	stream.pword(0) = stream.rdbuf();
	stream.register_callback(stream_cb, 0);
}

int
g3_istream_handle(std::istream &stream)
{
	std::streambuf* sbuf = stream.rdbuf();
	if (!sbuf)
		return -1;
	RemoteInputStreamBuffer* rbuf = dynamic_cast<RemoteInputStreamBuffer*>(sbuf);
	if (!rbuf)
		return -1;
	return rbuf->fd();
}

void
g3_ostream_to_path(std::ostream &stream, const std::string &path, bool append,
    size_t buffersize, const std::string &ext)
{
	reset_stream(stream);

	Codec codec = check_output_path(path, ext);
	if (append && codec != NONE)
		log_fatal("Cannot append to compressed file.");

	switch(codec) {
#ifdef ZLIB_FOUND
	case GZ:
		stream.rdbuf(new GZipEncoder(path, buffersize));
		break;
#endif
#ifdef BZIP2_FOUND
	case BZIP2:
		stream.rdbuf(new BZip2Encoder(path, buffersize));
		break;
#endif
#ifdef LZMA_FOUND
	case LZMA:
		stream.rdbuf(new LZMAEncoder(path, buffersize));
		break;
#endif
	default:
		stream.rdbuf(new OutputFileStreamCounter(path, buffersize, append));
		break;
	}

	stream.pword(0) = stream.rdbuf();
	stream.register_callback(stream_cb, 1);
}
