#include <G3Logging.h>
#include <dataio.h>
#include "compression.h"

#include <filesystem>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <streambuf>


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

class RemoteInputStreamBuffer : public std::streambuf {
public:
	RemoteInputStreamBuffer(const std::string &path, float timeout, size_t size)
	  : buffer_(size), bytes_(0) {
		fd_ = connect_remote(path, timeout);
		setg(buffer_.data(), buffer_.data(), buffer_.data());
	}

	~RemoteInputStreamBuffer() {
		close(fd_);
	}

	int fd() { return fd_; }

protected:
	int_type underflow() {
		if (gptr() < egptr())
			return traits_type::to_int_type(*gptr());

		ssize_t n = read(fd_, buffer_.data(), buffer_.size());
		if (n <= 0)
			return traits_type::eof();
		setg(buffer_.data(), buffer_.data(), buffer_.data() + n);
		return traits_type::to_int_type(*gptr());
	}

	std::streamsize xsgetn(char* s, std::streamsize n) {
		std::streamsize n_read = 0;
		while (n_read < n) {
			if (gptr() == egptr()) {
				if (underflow() == traits_type::eof())
					break;
			}

			std::streamsize remaining = n - n_read;
			std::streamsize available = egptr() - gptr();
			std::streamsize to_read = std::min(remaining, available);

			std::memcpy(s + n_read, gptr(), to_read);
			gbump(to_read);
			n_read += to_read;
		}
		bytes_ += n_read;
		return n_read;
	}

	std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,
	    std::ios_base::openmode mode) {
		// short-circuit for tellg
		if ((mode & std::ios_base::in) && off == 0 && way == std::ios_base::cur)
			return bytes_;
		log_fatal("Seek not implemented for remote stream");
	}

	std::streampos seekpos(std::streampos pos, std::ios_base::openmode mode) {
		log_fatal("Seek not implemented for remote stream");
	}

private:
	int fd_;
	std::vector<char> buffer_;
	size_t bytes_;
};


enum Codec {
  NONE = 0,
  GZ = 1,
  BZIP2 = 2,
  LZMA = 3,
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

	log_fatal("Invalid filename %s", path.c_str());
}

static Codec
check_input_path(const std::string &path, const std::string &ext)
{
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
	}

	return get_codec(path, ext);
}

class InputFileStreamCounter : public std::streambuf {
public:
	InputFileStreamCounter(const std::string& path, size_t size)
	  : buffer_(size), bytes_(0) {
		file_.open(path, std::ios::binary);
		if (!file_.is_open())
			log_fatal("Error opening file %s", path.c_str());
		file_.rdbuf()->pubsetbuf(buffer_.data(), buffer_.size());
	}

protected:
	int_type underflow() {
		return file_.rdbuf()->sgetc();
	}

	std::streamsize xsgetn(char* s, std::streamsize n) {
		std::streamsize nget = file_.rdbuf()->sgetn(s, n);
		if (nget > 0)
			bytes_ += nget;
		return nget;
	}

	std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,
	    std::ios_base::openmode mode) {
		if (!(mode & std::ios_base::in))
			log_fatal("Seek not implemented for output stream");
		// short-circuit for tellg
		if (off == 0 && way == std::ios_base::cur)
			return bytes_;
		std::streampos n = file_.rdbuf()->pubseekoff(off, way, mode);;
		if (n != std::streampos(std::streamoff(-1)))
			bytes_ = n;
		return n;
	}

	std::streampos seekpos(std::streampos pos, std::ios_base::openmode mode) {
		if (!(mode & std::ios_base::in))
			log_fatal("Seek not implemented for output stream");
		std::streampos n = file_.rdbuf()->pubseekpos(pos, mode);;
		if (n != std::streampos(std::streamoff(-1)))
			bytes_ = n;
		return n;
	}

private:
	std::ifstream file_;
	std::vector<char> buffer_;
	size_t bytes_;
};

class OutputFileStreamCounter : public std::streambuf {
public:
	OutputFileStreamCounter(const std::string& path, size_t size, bool append)
	  : buffer_(size), bytes_(0) {
		std::ios_base::openmode mode = std::ios::binary;
		if (append)
			mode |= std::ios::app;
		file_.open(path, mode);
		if (!file_.is_open())
			log_fatal("Error opening file %s", path.c_str());
		file_.rdbuf()->pubsetbuf(buffer_.data(), buffer_.size());
	}

protected:
	int_type overflow(int_type c) {
		if (file_.rdbuf()->sputc(c) != traits_type::eof()) {
			bytes_++;
			return c;
		}
		return traits_type::eof();
	}

	std::streamsize xsputn(const char* s, std::streamsize n) {
		std::streamsize nput = file_.rdbuf()->sputn(s, n);
		if (nput > 0)
			bytes_ += nput;
		return nput;
	}

	int sync() {
		return file_.rdbuf()->pubsync();
	}

	std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,
	    std::ios_base::openmode mode) {
		// short-circuit for tellp
		if ((mode & std::ios_base::out) && off == 0 && way == std::ios_base::cur)
			return bytes_;
		log_fatal("Seek not implemented for output stream");
	}

	std::streampos seekpos(std::streampos pos, std::ios_base::openmode mode) {
		log_fatal("Seek not implemented for output stream");
	}

private:
	std::ofstream file_;
	std::vector<char> buffer_;
	size_t bytes_;
};


static void
reset_stream(std::ios &stream)
{
	std::streambuf* sbuf = stream.rdbuf();
	if (sbuf) {
		sbuf->pubsync();
		delete sbuf;
	}
	stream.rdbuf(nullptr);
}

static void
stream_cb(std::ios::event ev, std::ios_base& stream, int index)
{
	std::ios* str = nullptr;

	switch (ev) {
	case std::ios::event::erase_event:
		str = dynamic_cast<std::ios*>(&stream);
		if (str)
			reset_stream(*str);
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

	// Figure out what kind of ultimate data source this is
	if (path.find("tcp://") == 0) {
		stream.rdbuf(new RemoteInputStreamBuffer(path, timeout, buffersize));
		return;
	}

	// Simple file case
	switch(check_input_path(path, ext)) {
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

	stream.register_callback(stream_cb, 1);
}
