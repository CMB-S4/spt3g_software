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

template<typename T, typename C>
class Decoder : public std::streambuf {
public:
	Decoder(std::istream& file, size_t size)
	  : file_(file), inbuf_(size), outbuf_(size) {
		setg(outbuf_.data(), outbuf_.data(), outbuf_.data());
	}

protected:
	int_type underflow() {
		if (gptr() < egptr())
			return traits_type::to_int_type(*gptr());

		stream_.avail_in = file_.read(inbuf_.data(), inbuf_.size()).gcount();
		if (stream_.avail_in <= 0)
			return traits_type::eof();
		stream_.next_in = reinterpret_cast<C*>(inbuf_.data());

		stream_.avail_out = outbuf_.size();
		stream_.next_out = reinterpret_cast<C*>(outbuf_.data());

		if (decode())
			return traits_type::eof();

		setg(outbuf_.data(), outbuf_.data(),
		    outbuf_.data() + outbuf_.size() - stream_.avail_out);
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

	std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,
	    std::ios_base::openmode mode) {
		log_fatal("Seek not implemented for compressed stream");
	}

	std::streampos seekpos(std::streampos pos, std::ios_base::openmode mode) {
		log_fatal("Seek not implemented for compressed stream");
	}

	virtual int decode() = 0;

	std::istream& file_;
	std::vector<char> inbuf_;
	std::vector<char> outbuf_;
	T stream_;
};

#ifdef ZLIB_FOUND
class GZipDecoder : public Decoder<z_stream, unsigned char> {
public:
	GZipDecoder(std::istream& file, size_t size) : Decoder(file, size)
	{
		stream_.zalloc = Z_NULL;
		stream_.zfree = Z_NULL;
		stream_.opaque = Z_NULL;
		stream_.avail_in = 0;
		stream_.next_in = Z_NULL;
		if (inflateInit2(&stream_, 16 + MAX_WBITS) != Z_OK)
			log_fatal("Error initializing gzip decoder: %s", stream_.msg);
	}

	~GZipDecoder() {
		inflateEnd(&stream_);
	}

protected:
	int decode() {
		int ret = inflate(&stream_, Z_NO_FLUSH);
		if (ret != Z_OK && ret != Z_STREAM_END) {
			log_error("Error running gzip decoder: %s", stream_.msg);
			return ret;
		}
		return 0;
	}
};
#endif

#ifdef BZIP2_FOUND
class BZip2Decoder : public Decoder<bz_stream, char> {
public:
	BZip2Decoder(std::istream& file, size_t size) : Decoder(file, size)
	{
		stream_.bzalloc = nullptr;
		stream_.bzfree = nullptr;
		stream_.opaque = nullptr;
		stream_.avail_in = 0;
		stream_.next_in = nullptr;
		if (BZ2_bzDecompressInit(&stream_, 0, 0) != BZ_OK)
			log_fatal("Error initializing bzip2 decoder");
	}

	~BZip2Decoder() {
		BZ2_bzDecompressEnd(&stream_);
	}

protected:
	int decode() {
		int ret = BZ2_bzDecompress(&stream_);
		if (ret != BZ_OK && ret != BZ_STREAM_END) {
			log_error("Error running bzip2 decoder");
			return ret;
		}
		return 0;
	}
};
#endif

#ifdef LZMA_FOUND
class LZMADecoder : public Decoder<lzma_stream, uint8_t> {
public:
	LZMADecoder(std::istream& file, size_t size) : Decoder(file, size)
        {
		stream_ = LZMA_STREAM_INIT;
		lzma_ret ret = lzma_stream_decoder(&stream_, UINT64_MAX,
		     LZMA_CONCATENATED);
		if (ret != LZMA_OK)
			log_fatal("Error initializing LZMA decoder.");
	}

	~LZMADecoder() {
		lzma_end(&stream_);
	}

protected:
	int decode() {
		lzma_ret ret = lzma_code(&stream_, LZMA_RUN);
		if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
			log_error("Error running LZMA decoder");
			return ret;
		}
		return 0;
	}
};
#endif

class InputStreamCounter : public std::streambuf {
public:
	InputStreamCounter(std::streambuf* buffer)
		: buffer_(buffer), bytes_(0) {
		if (!buffer_)
			log_fatal("Input file stream buffer required");
	}

protected:
	int_type underflow() {
		return buffer_->sgetc();
	}

	std::streamsize xsgetn(char* s, std::streamsize n) {
		std::streamsize nget = buffer_->sgetn(s, n);
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
		std::streampos n = buffer_->pubseekoff(off, way, mode);
		if (n != std::streampos(std::streamoff(-1)))
			bytes_ = n;
		return n;
	}

	std::streampos seekpos(std::streampos pos, std::ios_base::openmode mode) {
		if (!(mode & std::ios_base::in))
			log_fatal("Seek not implemented for output stream");
		std::streampos n = buffer_->pubseekpos(pos, mode);
		if (n != std::streampos(std::streamoff(-1)))
			bytes_ = n;
		return n;
	}

private:
	std::streambuf* buffer_;
	size_t bytes_;
};

void
g3_istream_from_path(std::istream &stream, const std::string &path, float timeout,
    size_t buffersize, const std::string &ext)
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

		switch(get_codec(path, ext)) {
#ifdef ZLIB_FOUND
		case GZ:
			sbuf = new GZipDecoder(*file, buffersize);
			break;
#endif
#ifdef BZIP2_FOUND
		case BZIP2:
			sbuf = new BZip2Decoder(*file, buffersize);
			break;
#endif
#ifdef LZMA_FOUND
		case LZMA:
			sbuf = new LZMADecoder(*file, buffersize);
			break;
#endif
		default:
			// Read buffer
			fbuf = new std::vector<char>(buffersize);
			file->rdbuf()->pubsetbuf(fbuf->data(), buffersize);
			sbuf = new InputStreamCounter(file->rdbuf());
			break;
		}
	}

	stream.rdbuf(sbuf);
	stream.pword(0) = file;
	stream.pword(1) = fbuf;
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
g3_istream_close(std::istream &stream)
{
	std::vector<char>* fbuf = static_cast<std::vector<char>*>(stream.pword(1));
	if (fbuf)
		delete fbuf;
	stream.pword(1) = nullptr;

	std::ifstream* file = static_cast<std::ifstream*>(stream.pword(0));
	if (file)
		delete file;
	stream.pword(0) = nullptr;

	std::streambuf* sbuf = stream.rdbuf();
	if (sbuf)
		delete sbuf;
	stream.rdbuf(nullptr);
}

class OutputStreamCounter : public std::streambuf {
public:
	OutputStreamCounter(std::streambuf* buffer)
	  : buffer_(buffer), bytes_(0) {
		if (!buffer_)
			log_fatal("Input file stream buffer required");
	}

protected:
	int_type overflow(int_type c) {
		if (buffer_->sputc(c) != traits_type::eof()) {
			bytes_++;
			return c;
		}
		return traits_type::eof();
	}

	std::streamsize xsputn(const char* s, std::streamsize n) {
		std::streamsize nput = buffer_->sputn(s, n);
		if (nput > 0)
			bytes_ += nput;
		return nput;
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
	std::streambuf* buffer_;
	size_t bytes_;
};

template <typename T, typename C>
class Encoder : public std::streambuf {
public:
	Encoder(std::ostream &file, size_t size)
	  : file_(file), inbuf_(size), outbuf_(size), bytes_(0) {}

protected:
	int_type overflow(int_type c) {
		if (pptr() && pbase()) {
			std::streamsize len = pptr() - pbase();
			stream_.avail_in = len;
			stream_.next_in = reinterpret_cast<C*>(pbase());
			encode_all();
		}
		if (c != traits_type::eof()) {
			inbuf_[0] = traits_type::to_char_type(c);
			stream_.avail_in = 1;
			stream_.next_in = reinterpret_cast<C*>(inbuf_.data());
			encode_all();
		}
		setp(inbuf_.data(), inbuf_.data() + inbuf_.size());
		return traits_type::not_eof(c);
	}

	std::streamsize xsputn(const char* s, std::streamsize n) {
		stream_.avail_in = n;
		stream_.next_in = reinterpret_cast<C*>(const_cast<char*>(s));
		encode_all();
		return n;
	}

	int sync() {
		stream_.avail_in = 0;
		encode_all(true);
		file_.flush();
		return 0;
	}

	void encode_all(bool flush=false) {
		do {
			stream_.avail_out = outbuf_.size();
			stream_.next_out = reinterpret_cast<C*>(outbuf_.data());
			if (encode(flush))
				return;
			size_t n =  outbuf_.size() - stream_.avail_out;
			bytes_ += n;
			file_.write(outbuf_.data(), n);
		} while (stream_.avail_out == 0);
	}

	virtual int encode(bool flush) = 0;

	std::streampos seekoff(std::streamoff off, std::ios_base::seekdir way,
	    std::ios_base::openmode mode) {
		// short-circuit for tellp
		if ((mode & std::ios_base::out) && off == 0 && way == std::ios_base::cur)
			return bytes_;
		log_fatal("Seek not implemented for compressed stream");
	}

	std::streampos seekpos(std::streampos pos, std::ios_base::openmode mode) {
		log_fatal("Seek not implemented for compressed stream");
	}

	std::ostream &file_;
	std::vector<char> inbuf_;
	std::vector<char> outbuf_;
	size_t bytes_;
	T stream_;
};

#ifdef ZLIB_FOUND
class GZipEncoder : public Encoder<z_stream, unsigned char> {
public:
	GZipEncoder(std::ostream &file, size_t size) : Encoder(file, size)
	{
		stream_.zalloc = Z_NULL;
		stream_.zfree = Z_NULL;
		stream_.opaque = Z_NULL;
		if (deflateInit2(&stream_, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
		    16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
			log_fatal("Error initializing gzip encoder: %s", stream_.msg);
	}

	~GZipEncoder() {
		deflateEnd(&stream_);
	}

protected:
	int encode(bool flush) {
		int ret = deflate(&stream_, flush ? Z_FINISH : Z_NO_FLUSH);
		if (ret == Z_STREAM_ERROR) {
			log_error("Error running gzip encoder: %s", stream_.msg);
			return ret;
		}
		return 0;
	}
};
#endif

#ifdef BZIP2_FOUND
class BZip2Encoder : public Encoder<bz_stream, char> {
public:
	BZip2Encoder(std::ostream &file, size_t size) : Encoder(file, size)
	{
		stream_.bzalloc = nullptr;
		stream_.bzfree = nullptr;
		stream_.opaque = nullptr;
		if (BZ2_bzCompressInit(&stream_, 9, 0, 0) != BZ_OK)
			log_fatal("Error initializing bzip2 encoder");
	}

	~BZip2Encoder() {
		BZ2_bzCompressEnd(&stream_);
	}

protected:
	int encode(bool flush) {
		int ret = BZ2_bzCompress(&stream_, flush ? BZ_FINISH : BZ_RUN);
		if (ret == BZ_SEQUENCE_ERROR) {
			log_error("Error running bzip2 encoder");
			return ret;
		}
		return 0;
	}
};
#endif

#ifdef LZMA_FOUND
class LZMAEncoder : public Encoder<lzma_stream, uint8_t> {
public:
	LZMAEncoder(std::ostream &file, size_t size) : Encoder(file, size)
	{
		stream_ = LZMA_STREAM_INIT;
		lzma_ret ret = lzma_easy_encoder(&stream_, LZMA_PRESET_DEFAULT,
		    LZMA_CHECK_CRC64);
		if (ret != LZMA_OK)
			log_fatal("Error initializing LZMA encoder.");
	}

	~LZMAEncoder() {
		lzma_end(&stream_);
	}

protected:
	int encode(bool flush) {
		lzma_action action = flush ? LZMA_FINISH : LZMA_RUN;
		lzma_ret ret = lzma_code(&stream_, action);
		if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
			log_error("Error running LZMA encoder");
			return ret;
		}
		return 0;
	}
};
#endif

void
g3_ostream_to_path(std::ostream &stream, const std::string &path, bool append,
    size_t buffersize, const std::string &ext)
{
	Codec codec = get_codec(path, ext);
	std::ios_base::openmode mode = std::ios::binary;
	if (append) {
		if (codec != NONE)
			log_fatal("Cannot append to compressed file.");
		mode |= std::ios::app;
	}
	std::ofstream *file = new std::ofstream(path, mode);
	if (!file->is_open()) {
		delete file;
		log_fatal("Could not open file %s", path.c_str());
	}

	std::streambuf *sbuf = nullptr;

	switch(codec) {
#ifdef ZLIB_FOUND
	case GZ:
		sbuf = new GZipEncoder(*file, buffersize);
		break;
#endif
#ifdef BZIP2_FOUND
	case BZIP2:
		sbuf = new BZip2Encoder(*file, buffersize);
		break;
#endif
#ifdef LZMA_FOUND
	case LZMA:
		sbuf = new LZMAEncoder(*file, buffersize);
		break;
#endif
	default:
		sbuf = new OutputStreamCounter(file->rdbuf());
		break;
	}

	stream.rdbuf(sbuf);
	stream.pword(0) = file;
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
	if (file)
		delete file;
	stream.pword(0) = nullptr;

	std::streambuf* sbuf = stream.rdbuf();
	if (sbuf)
		delete sbuf;
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
