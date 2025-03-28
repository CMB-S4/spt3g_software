#ifndef _G3_COMPRESSION_H
#define _G3_COMPRESSION_H

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

#include <fstream>
#include <streambuf>

template<typename T, typename C>
class Decoder : public std::streambuf {
public:
	Decoder(const std::string& path, size_t size)
	  : inbuf_(size), outbuf_(size) {
		file_.open(path, std::ios::binary);
		if (!file_.is_open())
			log_fatal("Could not open file %s", path.c_str());
		setg(outbuf_.data(), outbuf_.data(), outbuf_.data());
	}

	virtual ~Decoder() {
		if (file_.is_open())
			file_.close();
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

	std::ifstream file_;
	std::vector<char> inbuf_;
	std::vector<char> outbuf_;
	T stream_;
};

#ifdef ZLIB_FOUND
class GZipDecoder : public Decoder<z_stream, unsigned char> {
public:
	GZipDecoder(const std::string& path, size_t size);
	~GZipDecoder();

protected:
	int decode();
};
#endif

#ifdef BZIP2_FOUND
class BZip2Decoder : public Decoder<bz_stream, char> {
public:
	BZip2Decoder(const std::string& path, size_t size);
	~BZip2Decoder();

protected:
	int decode();
};
#endif

#ifdef LZMA_FOUND
class LZMADecoder : public Decoder<lzma_stream, uint8_t> {
public:
	LZMADecoder(const std::string& path, size_t size);
	~LZMADecoder();

protected:
	int decode();
};
#endif

template <typename T, typename C>
class Encoder : public std::streambuf {
public:
	Encoder(const std::string& path, size_t size)
	  : inbuf_(size), outbuf_(size), bytes_(0) {
		file_.open(path, std::ios::binary);
		if (!file_.is_open())
			log_fatal("Could not open file %s", path.c_str());
	}

	virtual ~Encoder() {
		if (file_.is_open()) {
			file_.flush();
			file_.close();
		}
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
		return file_.rdbuf()->pubsync();
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

	std::ofstream file_;
	std::vector<char> inbuf_;
	std::vector<char> outbuf_;
	size_t bytes_;
	T stream_;
};

#ifdef ZLIB_FOUND
class GZipEncoder : public Encoder<z_stream, unsigned char> {
public:
	GZipEncoder(const std::string& path, size_t size);
	~GZipEncoder();

protected:
	int encode(bool flush);
};
#endif

#ifdef BZIP2_FOUND
class BZip2Encoder : public Encoder<bz_stream, char> {
public:
	BZip2Encoder(const std::string& path, size_t size);
	~BZip2Encoder();

protected:
	int encode(bool flush);
};
#endif

#ifdef LZMA_FOUND
class LZMAEncoder : public Encoder<lzma_stream, uint8_t> {
public:
	LZMAEncoder(const std::string& path, size_t size);
	~LZMAEncoder();

protected:
	int encode(bool flush);
};
#endif

#endif
