#include <G3Logging.h>
#include <dataio.h>

#ifdef ZLIB_FOUND
#include <zlib.h>
#endif
#ifdef BZIP2_FOUND
#include <bzlib.h>
#endif
#ifdef LZMA_FOUND
#include <lzma.h>
#endif

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
connect_remote(const std::string &path, int timeout)
{
	// TCP Socket. Two syntaxes:
	// - tcp://host:port -> connect to "host" on "port" and read
	//   until EOF
	// - tcp://*:port -> listen on "port" for the first connection
	//   and read until EOF

	std::string host = path.substr(path.find("://") + 3);
	if (host.find(":") == host.npos)
		log_fatal("Could not open URL %s: unspecified port",
		    path.c_str());
	std::string port = host.substr(host.find(":") + 1);
	host = host.substr(0, host.find(":"));

	log_debug("Opening connection to %s, port %s", host.c_str(),
	    port.c_str());

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
	} else {
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
	}

	return fd;
}

class RemoteInputStreamBuffer : public std::streambuf {
public:
	RemoteInputStreamBuffer(const std::string &path, int timeout, size_t size)
	  : buffer_(size) {
		fd_ = connect_remote(path, timeout);
		setg(buffer_.data(), buffer_.data(), buffer_.data());
	}

	~RemoteInputStreamBuffer() {
		close(fd_);
	}

	int fd() { return fd_; }

protected:
	int_type underflow() {
		ssize_t n = read(fd_, buffer_.data(), buffer_.size());
		if (n <= 0)
			return traits_type::eof();
		setg(buffer_.data(), buffer_.data(), buffer_.data() + n);
		return traits_type::to_int_type(buffer_[0]);
	}

private:
	int fd_;
	std::vector<char> buffer_;
};

class GZipDecompressor : public std::streambuf {
public:
	GZipDecompressor(std::istream& file, size_t size)
#ifdef ZLIB_FOUND
	  : file_(file), buffer_(size) {
		(void) file_;
		log_fatal("Not implmented");
	}

private:
	std::istream& file_;
	std::vector<char> buffer_;
#else
	{
		log_fatal("Library not compiled with gzip support.");
	}
#endif
};

class BZip2Decompressor : public std::streambuf {
public:
	BZip2Decompressor(std::istream& file, size_t size)
#ifdef BZIP2_FOUND
	  : file_(file), buffer_(size) {
		(void) file_;
		log_fatal("Not implmented");
	}

private:
	std::istream& file_;
	std::vector<char> buffer_;
#else
	{
		log_fatal("Library not compiled with bzip2 support.");
	}
#endif
};

class LZMADecompressor : public std::streambuf {
public:
	LZMADecompressor(std::istream& file, size_t size)
#ifdef LZMA_FOUND
	  : file_(file), buffer_(size) {
		(void) file_;
		log_fatal("Not implmented");
	}

private:
	std::istream& file_;
	std::vector<char> buffer_;
#else
	{
		log_fatal("Library not compiled with LZMA support.");
	}
#endif
};

void
g3_istream_from_path(std::istream &stream, const std::string &path, float timeout,
    size_t buffersize)
{
	g3_istream_close(stream);

	std::streambuf *sbuf = nullptr;
	std::ifstream *file = nullptr;
	std::vector<char> *fbuf = nullptr;

	// Figure out what kind of ultimate data source this is
	if (path.find("tcp://") == 0) {
		sbuf = new RemoteInputStreamBuffer(path, timeout, buffersize);
	} else {
		// Simple file case
		file = new std::ifstream(path, std::ios::binary);
		if (!file->is_open()) {
			delete file;
			log_fatal("Could not open file %s", path.c_str());
		}

		// Read buffer
		fbuf = new std::vector<char>(buffersize);
		file->rdbuf()->pubsetbuf(fbuf->data(), buffersize);

		if (path.size() > 3 && !path.compare(path.size() - 3, 3, ".gz")) {
			sbuf = new GZipDecompressor(*file, buffersize);
		} else if (path.size() > 4 && !path.compare(path.size() - 4, 4, ".bz2")) {
			sbuf = new BZip2Decompressor(*file, buffersize);
		} else if (path.size() > 3 && !path.compare(path.size() - 3, 3, ".xz")) {
			sbuf = new LZMADecompressor(*file, buffersize);
		} else {
			sbuf = file->rdbuf();
		}
	}

	stream.rdbuf(sbuf);
	stream.pword(0) = file;
	stream.pword(1) = fbuf;
	stream.pword(2) = sbuf;
}

int
g3_istream_handle(std::istream &stream)
{
	std::streambuf* sbuf = static_cast<std::streambuf*>(stream.pword(2));
	if (!sbuf)
		return -1;
	RemoteInputStreamBuffer* rbuf = dynamic_cast<RemoteInputStreamBuffer*>(sbuf);
	if (!rbuf)
		return -1;
	return rbuf->fd();
}

void
g3_istream_close(std::istream &stream)
{
	std::ifstream* file = static_cast<std::ifstream*>(stream.pword(0));
	std::vector<char>* fbuf = static_cast<std::vector<char>*>(stream.pword(1));
	std::streambuf* sbuf = static_cast<std::streambuf*>(stream.pword(2));

	if (sbuf && (!file || sbuf != file->rdbuf()))
		delete sbuf;
	stream.pword(2) = nullptr;

	if (fbuf)
		delete fbuf;
	stream.pword(1) = nullptr;

	if (file)
		delete file;
	stream.pword(0) = nullptr;

	stream.rdbuf(nullptr);
}

class Compressor : public std::streambuf {
public:
	Compressor(std::ostream &file, size_t size)
	  : file_(file), buffer_(size), idx_(0) {}

protected:
	int_type overflow(int_type c) {
		if (c == traits_type::eof())
			return sync() == 0 ? 0 : traits_type::eof();

		buffer_[idx_] = traits_type::to_char_type(c);
		idx_++;
		if (idx_ == buffer_.size()) {
			if (sync())
				return traits_type::eof();
		}
		return c;
	}

	int sync() {
		if (idx_ > 0) {
			if (compress())
				return -1;
			idx_ = 0;
		}

		if (finalize())
			return -1;

		return 0;
	}

	virtual int compress() { return 0; };
	virtual int finalize() { return 0; };

	std::ostream &file_;
	std::vector<char> buffer_;
	size_t idx_;
};

class GZipCompressor : public Compressor {
public:
	GZipCompressor(std::ostream &file, size_t size)
	  : Compressor(file, size)
#ifdef ZLIB_FOUND
	{
		(void) file_;
		log_fatal("Not implemented");
	}
#else
	{
		log_fatal("Library not implemented with gzip support.");
	}
#endif
};

class BZip2Compressor : public Compressor {
public:
	BZip2Compressor(std::ostream &file, size_t size)
	  : Compressor(file, size)
#ifdef BZIP2_FOUND
	{
		(void) file_;
		log_fatal("Not implemented");
	}
#else
	{
		log_fatal("Library not implemented with bzip2 support.");
	}
#endif
};

class LZMACompressor : public Compressor {
public:
	LZMACompressor(std::ostream &file, size_t size)
	  : Compressor(file, size)
#ifdef LZMA_FOUND
	{
		(void) file_;
		log_fatal("Not implemented");
	}
#else
	{
		log_fatal("Library not implemented with LZMA support.");
	}
#endif
};

class OutputStreamCounter : public std::streambuf {
public:
	OutputStreamCounter(std::streambuf* buffer)
	  : buffer_(buffer), bytes_(0) {}

	std::streampos bytes() const { return bytes_; }

protected:
	int_type overflow(int_type c) {
		if (buffer_ && buffer_->sputc(c) != traits_type::eof()) {
			bytes_++;
			return c;
		}
		return traits_type::eof();
	}

	std::streamsize xsputn(const char* s, std::streamsize n) {
		if (!buffer_)
			return 0;

		std::streamsize nput = buffer_->sputn(s, n);
		if (nput > 0)
			bytes_ += nput;
		return nput;
	}

private:
	std::streambuf* buffer_;
	size_t bytes_;
};

void
g3_ostream_to_path(std::ostream &stream, const std::string &path, bool append,
    bool counter, size_t buffersize)
{
	std::ios_base::openmode mode = std::ios::binary;
	if (append)
		mode |= std::ios::app;
	std::ofstream *file = new std::ofstream(path, mode);
	if (!file->is_open()) {
		delete file;
		log_fatal("Could not open file %s", path.c_str());
	}

	std::streambuf *sbuf = nullptr;

	if (path.size() > 3 && !path.compare(path.size() - 3, 3, ".gz")) {
		if (append)
			log_fatal("Cannot append to compressed file.");
		sbuf = new GZipCompressor(*file, buffersize);
	} else if (path.size() > 4 && !path.compare(path.size() - 4, 4, ".bz2")) {
		if (append)
			log_fatal("Cannot append to compressed file.");
		sbuf = new BZip2Compressor(*file, buffersize);
	} else if (path.size() > 3 && !path.compare(path.size() - 3, 3, ".xz")) {
		if (append)
			log_fatal("Cannot append to compressed file.");
		sbuf = new LZMACompressor(*file, buffersize);
	} else {
		sbuf = file->rdbuf();
	}

	std::streambuf *cbuf = nullptr;
	if (counter) {
		cbuf = new OutputStreamCounter(sbuf);
		stream.rdbuf(cbuf);
	} else {
		stream.rdbuf(sbuf);
	}

	stream.pword(0) = file;
	stream.pword(1) = sbuf;
	stream.pword(2) = cbuf;
}

size_t
g3_ostream_count(std::ostream &stream)
{
	OutputStreamCounter* cbuf = static_cast<OutputStreamCounter*>(stream.pword(2));
	if (!cbuf)
		log_fatal("Could not get stream counter");

	return cbuf->bytes();
}

void
g3_ostream_flush(std::ostream &stream)
{
	stream.flush();

	std::streambuf* sbuf = stream.rdbuf();
	if (sbuf)
		sbuf->pubsync();

	std::ofstream* file = static_cast<std::ofstream*>(stream.pword(0));
	if (file)
		file->flush();
}

void
g3_ostream_close(std::ostream &stream)
{
	g3_ostream_flush(stream);

	std::ofstream* file = static_cast<std::ofstream*>(stream.pword(0));
	std::streambuf* sbuf = static_cast<std::streambuf*>(stream.pword(1));
	if (sbuf) {
		if (!file || (sbuf != file->rdbuf()))
			delete sbuf;
	}
	stream.pword(1) = nullptr;

	if (file)
		delete file;
	stream.pword(0) = nullptr;

	OutputStreamCounter* cbuf = static_cast<OutputStreamCounter*>(stream.pword(2));
	if (cbuf)
		delete cbuf;
	stream.pword(2) = nullptr;

	stream.rdbuf(nullptr);
}

void
g3_check_input_path(const std::string &path)
{
	if (path.find("://") != path.npos)
		return;

	std::filesystem::path fpath(path);
	if (!std::filesystem::exists(fpath) ||
	    !std::filesystem::is_regular_file(fpath))
		log_fatal("Could not find file %s", path.c_str());
}

void
g3_check_output_path(const std::string &path)
{
	std::filesystem::path fpath(path);

	if (fpath.empty())
		log_fatal("Empty file path");

	if (!fpath.has_parent_path())
		return;

	if (!std::filesystem::exists(fpath.parent_path()))
		log_fatal("Parent path does not exist: %s",
		    fpath.parent_path().string().c_str());
}
