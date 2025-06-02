#ifndef _G3_STREAMS_H
#define _G3_STREAMS_H

#include <G3Logging.h>

#ifdef ZLIB_FOUND
#include <zlib.h>
#endif
#ifdef BZIP2_FOUND
#include <bzlib.h>
#endif
#ifdef LZMA_FOUND
#include <lzma.h>
#endif

#include <unistd.h>
#include <fstream>
#include <streambuf>

class RemoteInputStreamBuffer : public std::streambuf {
public:
	RemoteInputStreamBuffer(int fd, size_t size)
	  : fd_(fd), buffer_(new char[size]), bsize_(size), bytes_(0) {
		setg(&buffer_[0], &buffer_[0], &buffer_[0]);
	}

	~RemoteInputStreamBuffer() {
		close(fd_);
	}

	int fd() { return fd_; }

protected:
	int_type underflow() {
		if (gptr() < egptr())
			return traits_type::to_int_type(*gptr());

		ssize_t n = read(fd_, &buffer_[0], bsize_);
		if (n <= 0)
			return traits_type::eof();
		setg(&buffer_[0], &buffer_[0], &buffer_[0] + n);
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
	std::unique_ptr<char[]> buffer_;
	size_t bsize_;
	size_t bytes_;
};

class InputFileStreamCounter : public std::filebuf {
public:
	InputFileStreamCounter(const std::string& path, size_t size)
	  : buffer_(new char[size]), bytes_(0) {
		open(path, std::ios::in | std::ios::binary);
		if (!is_open())
			log_fatal("Error opening file %s", path.c_str());
		pubsetbuf(&buffer_[0], size);
	}

protected:
	std::streamsize xsgetn(char* s, std::streamsize n) {
		std::streamsize nget = std::filebuf::xsgetn(s, n);
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
		bytes_ = std::filebuf::seekoff(off, way, mode);
		return bytes_;
	}

	std::streampos seekpos(std::streampos pos, std::ios_base::openmode mode) {
		bytes_ = std::filebuf::seekpos(pos, mode);
		return bytes_;
	}

private:
	std::unique_ptr<char []> buffer_;
	size_t bytes_;
};

class OutputFileStreamCounter : public std::filebuf {
public:
	OutputFileStreamCounter(const std::string& path, size_t size, bool append)
	  : buffer_(new char[size]), bytes_(0) {
		std::ios_base::openmode mode = std::ios::out | std::ios::binary;
		if (append)
			mode |= std::ios::app;
		open(path, mode);
		if (!is_open())
			log_fatal("Error opening file %s", path.c_str());
		// Update bytes counter so that tellp works correctly
		if (append)
			bytes_ = std::filebuf::seekoff(0, std::ios_base::cur,
			    std::ios::out);
		pubsetbuf(&buffer_[0], size);
	}

protected:
	int_type overflow(int_type c) {
		int_type r = std::filebuf::overflow(c);
		if (r != traits_type::eof())
			bytes_++;
		return r;
	}

	std::streamsize xsputn(const char* s, std::streamsize n) {
		std::streamsize nput = std::filebuf::xsputn(s, n);
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
	std::unique_ptr<char []> buffer_;
	size_t bytes_;
};

template<typename T, typename C>
class Decoder : public std::streambuf {
public:
	Decoder(const std::string& path, size_t size)
	  : inbuf_(new char[size]), outbuf_(new char[size]), bsize_(size) {
		file_.open(path, std::ios::binary);
		if (!file_.is_open())
			log_fatal("Could not open file %s", path.c_str());
		setg(&outbuf_[0], &outbuf_[0], &outbuf_[0]);
	}

protected:
	int_type underflow() {
		if (gptr() < egptr())
			return traits_type::to_int_type(*gptr());

		if (stream_.avail_in == 0) {
			if (file_.eof())
				return traits_type::eof();
			stream_.avail_in = file_.read(&inbuf_[0], bsize_).gcount();
			if (stream_.avail_in <= 0)
				return traits_type::eof();
			stream_.next_in = reinterpret_cast<C*>(&inbuf_[0]);
		}

		stream_.avail_out = bsize_;
		stream_.next_out = reinterpret_cast<C*>(&outbuf_[0]);

		if (decode())
			return traits_type::eof();
		if (stream_.avail_out == bsize_)
			return traits_type::eof();

		setg(&outbuf_[0], &outbuf_[0],
		    &outbuf_[0] + bsize_ - stream_.avail_out);
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

	std::ifstream file_;
	std::unique_ptr<char []> inbuf_;
	std::unique_ptr<char []> outbuf_;
	size_t bsize_;
	T stream_;
};

template <typename T, typename C>
class Encoder : public std::streambuf {
public:
	Encoder(const std::string& path, size_t size)
	  : inbuf_(new char[size]), outbuf_(new char[size]), bsize_(size), bytes_(0) {
		file_.open(path, std::ios::binary);
		if (!file_.is_open())
			log_fatal("Could not open file %s", path.c_str());
	}

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
			stream_.next_in = reinterpret_cast<C*>(&inbuf_[0]);
			encode_all();
		}
		setp(&inbuf_[0], &inbuf_[0] + bsize_);
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
		return file_.rdbuf()->pubsync();
	}

	void encode_all(bool flush=false) {
		do {
			stream_.avail_out = bsize_;
			stream_.next_out = reinterpret_cast<C*>(&outbuf_[0]);
			if (encode(flush))
				return;
			size_t n =  bsize_ - stream_.avail_out;
			bytes_ += n;
			file_.write(&outbuf_[0], n);
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

	std::ofstream file_;
	std::unique_ptr<char []> inbuf_;
	std::unique_ptr<char []> outbuf_;
	size_t bsize_;
	size_t bytes_;
	T stream_;
};

#define CODEC(name, stype, ctype) \
class name##Decoder : public Decoder<stype, ctype> { \
public: \
	name##Decoder(const std::string& path, size_t size); \
	~name##Decoder(); \
protected: \
	int decode(); \
}; \
class name##Encoder : public Encoder<stype, ctype> { \
public: \
	name##Encoder(const std::string& path, size_t size); \
	~name##Encoder(); \
protected: \
	int encode(bool flush); \
}

#ifdef ZLIB_FOUND
CODEC(GZip, z_stream, unsigned char);
#endif
#ifdef BZIP2_FOUND
CODEC(BZip2, bz_stream, char);
#endif
#ifdef LZMA_FOUND
CODEC(LZMA, lzma_stream, uint8_t);
#endif

#endif
