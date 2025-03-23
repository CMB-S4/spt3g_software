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
		return n_read;
	}

private:
	int fd_;
	std::vector<char> buffer_;
};

class Decompressor : public std::streambuf {
public:
	Decompressor(std::istream& file, size_t size)
	  : file_(file), inbuf_(size), outbuf_(size) {
		setg(outbuf_.data(), outbuf_.data(), outbuf_.data());
	}

protected:
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
		return n_read;
	}

	std::istream& file_;
	std::vector<char> inbuf_;
	std::vector<char> outbuf_;
};

class GZipDecompressor : public Decompressor {
public:
	GZipDecompressor(std::istream& file, size_t size)
	  : Decompressor(file, size)
#ifdef ZLIB_FOUND
	{
		stream_.zalloc = Z_NULL;
		stream_.zfree = Z_NULL;
		stream_.opaque = Z_NULL;
		stream_.avail_in = 0;
		stream_.next_in = Z_NULL;
		inflateInit2(&stream_, 16 + MAX_WBITS);
	}

	~GZipDecompressor() {
		inflateEnd(&stream_);
	}

protected:
	int_type underflow() {
		if (gptr() < egptr())
			return traits_type::to_int_type(*gptr());

		stream_.avail_in = file_.read(inbuf_.data(), inbuf_.size()).gcount();
		if (stream_.avail_in <= 0)
			return traits_type::eof();
		stream_.next_in = reinterpret_cast<unsigned char*>(inbuf_.data());

		stream_.avail_out = outbuf_.size();
		stream_.next_out = reinterpret_cast<unsigned char*>(outbuf_.data());

		int ret = inflate(&stream_, Z_NO_FLUSH);
		if (ret != Z_OK && ret != Z_STREAM_END)
			return traits_type::eof();

		setg(outbuf_.data(), outbuf_.data(),
		    outbuf_.data() + outbuf_.size() - stream_.avail_out);
		return traits_type::to_int_type(*gptr());
	}

private:
	z_stream stream_;
#else
	{
		log_fatal("Library not compiled with gzip support.");
	}
#endif
};

class BZip2Decompressor : public Decompressor {
public:
	BZip2Decompressor(std::istream& file, size_t size)
	  : Decompressor(file, size)
#ifdef BZIP2_FOUND
	{
		stream_.bzalloc = nullptr;
		stream_.bzfree = nullptr;
		stream_.opaque = nullptr;
		stream_.avail_in = 0;
		stream_.next_in = nullptr;
		BZ2_bzDecompressInit(&stream_, 0, 0);
	}

	~BZip2Decompressor() {
		BZ2_bzDecompressEnd(&stream_);
	}

protected:
	int underflow() override {
		if (gptr() < egptr())
			return traits_type::to_int_type(*gptr());

		stream_.avail_in = file_.read(inbuf_.data(), inbuf_.size()).gcount();
		if (stream_.avail_in <= 0)
			return traits_type::eof();
		stream_.next_in = inbuf_.data();

		stream_.avail_out = outbuf_.size();
		stream_.next_out = outbuf_.data();

		int ret = BZ2_bzDecompress(&stream_);
		if (ret != BZ_OK && ret != BZ_STREAM_END)
			return traits_type::eof();

		setg(outbuf_.data(), outbuf_.data(),
		    outbuf_.data() + outbuf_.size() - stream_.avail_out);
		return traits_type::to_int_type(*gptr());
	}

private:
	bz_stream stream_;
#else
	{
		log_fatal("Library not compiled with bzip2 support.");
	}
#endif
};

class LZMADecompressor : public Decompressor {
public:
	LZMADecompressor(std::istream& file, size_t size)
	  : Decompressor(file, size)
#ifdef LZMA_FOUND
        {
		(void) file_;
		log_fatal("Not implmented");
	}
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

		if (path.size() > 3 && !path.compare(path.size() - 3, 3, ".gz")) {
			sbuf = new GZipDecompressor(*file, buffersize);
		} else if (path.size() > 4 && !path.compare(path.size() - 4, 4, ".bz2")) {
			sbuf = new BZip2Decompressor(*file, buffersize);
		} else if (path.size() > 3 && !path.compare(path.size() - 3, 3, ".xz")) {
			sbuf = new LZMADecompressor(*file, buffersize);
		} else {
			// Read buffer
			fbuf = new std::vector<char>(buffersize);
			file->rdbuf()->pubsetbuf(fbuf->data(), buffersize);
			sbuf = file->rdbuf();
		}
	}

	stream.rdbuf(sbuf);
	stream.pword(0) = file;
	stream.pword(1) = sbuf;
	stream.pword(2) = fbuf;
}

int
g3_istream_handle(std::istream &stream)
{
	std::streambuf* sbuf = static_cast<std::streambuf*>(stream.pword(1));
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
	std::streambuf* sbuf = static_cast<std::streambuf*>(stream.pword(1));
	std::vector<char>* fbuf = static_cast<std::vector<char>*>(stream.pword(2));

	if (fbuf)
		delete fbuf;
	stream.pword(2) = nullptr;

	if (sbuf && (!file || sbuf != file->rdbuf()))
		delete sbuf;
	stream.pword(1) = nullptr;

	if (file)
		delete file;
	stream.pword(0) = nullptr;

	stream.rdbuf(nullptr);
}

class Compressor : public std::streambuf {
public:
	Compressor(std::ostream &file, size_t size)
	  : file_(file), inbuf_(size), outbuf_(size) {}

protected:
	std::ostream &file_;
	std::vector<char> inbuf_;
	std::vector<char> outbuf_;
};

class GZipCompressor : public Compressor {
public:
	GZipCompressor(std::ostream &file, size_t size)
	  : Compressor(file, size)
#ifdef ZLIB_FOUND
	{
		stream_.zalloc = Z_NULL;
		stream_.zfree = Z_NULL;
		stream_.opaque = Z_NULL;
		deflateInit2(&stream_, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
		    16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY);
	}

	~GZipCompressor() {
		deflateEnd(&stream_);
	}

protected:
	int_type overflow(int_type c) {
		if (pptr() && pbase()) {
			std::streamsize len = pptr() - pbase();
			stream_.avail_in = len;
			stream_.next_in = reinterpret_cast<unsigned char*>(pbase());
			compress();
		}
		if (c != traits_type::eof()) {
			inbuf_[0] = traits_type::to_char_type(c);
			stream_.avail_in = 1;
			stream_.next_in = reinterpret_cast<unsigned char*>(inbuf_.data());
			compress();
		}
		setp(inbuf_.data(), inbuf_.data() + inbuf_.size());
		return traits_type::not_eof(c);
	}

	std::streamsize xsputn(const char* s, std::streamsize n) {
		stream_.avail_in = n;
		stream_.next_in = reinterpret_cast<unsigned char*>(const_cast<char*>(s));
		compress();
		return n;
	}

	int sync() {
		stream_.avail_in = 0;
		compress(Z_FINISH);
		file_.flush();
		return 0;
	}

private:
	void compress(int flush=Z_NO_FLUSH) {
		do {
			stream_.avail_out = outbuf_.size();
			stream_.next_out = reinterpret_cast<unsigned char*>(outbuf_.data());
			deflate(&stream_, flush);
			file_.write(outbuf_.data(), outbuf_.size() - stream_.avail_out);
		} while (stream_.avail_out == 0);
	}

	z_stream stream_;
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
		stream_.bzalloc = nullptr;
		stream_.bzfree = nullptr;
		stream_.opaque = nullptr;
		BZ2_bzCompressInit(&stream_, 9, 0, 0);
	}

	~BZip2Compressor() {
		BZ2_bzCompressEnd(&stream_);
	}

protected:
	int_type overflow(int_type c) {
		if (pptr() && pbase()) {
			std::streamsize len = pptr() - pbase();
			stream_.avail_in = len;
			stream_.next_in = pbase();
			compress();
		}
		if (c != traits_type::eof()) {
			inbuf_[0] = traits_type::to_char_type(c);
			stream_.avail_in = 1;
			stream_.next_in = inbuf_.data();
			compress();
		}
		setp(inbuf_.data(), inbuf_.data() + inbuf_.size());
		return traits_type::not_eof(c);
	}

	std::streamsize xsputn(const char* s, std::streamsize n) {
		stream_.avail_in = n;
		stream_.next_in = const_cast<char*>(s);
		compress();
		return n;
	}

	int sync() {
		stream_.avail_in = 0;
		compress(BZ_FINISH);
		file_.flush();
		return 0;
	}

private:
	void compress(int flush=BZ_RUN) {
		do {
			stream_.avail_out = outbuf_.size();
			stream_.next_out = outbuf_.data();
			BZ2_bzCompress(&stream_, flush);
			file_.write(outbuf_.data(), outbuf_.size() - stream_.avail_out);
		} while (stream_.avail_out == 0);
	}

	bz_stream stream_;
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
