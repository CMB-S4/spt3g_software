#include "compression.h"

#ifdef ZLIB_FOUND
GZipDecoder::GZipDecoder(const std::string& path, size_t size) : Decoder(path, size)
{
	stream_.zalloc = Z_NULL;
	stream_.zfree = Z_NULL;
	stream_.opaque = Z_NULL;
	stream_.avail_in = 0;
	stream_.next_in = Z_NULL;
	if (inflateInit2(&stream_, 16 + MAX_WBITS) != Z_OK)
		log_fatal("Error initializing gzip decoder: %s", stream_.msg);
}

GZipDecoder::~GZipDecoder()
{
	inflateEnd(&stream_);
}

int GZipDecoder::decode()
{
	int ret = inflate(&stream_, Z_NO_FLUSH);
	if (ret != Z_OK && ret != Z_STREAM_END) {
		log_error("Error running gzip decoder: %s", stream_.msg);
		return ret;
	}
	return 0;
}

GZipEncoder::GZipEncoder(const std::string& path, size_t size) : Encoder(path, size)
{
	stream_.zalloc = Z_NULL;
	stream_.zfree = Z_NULL;
	stream_.opaque = Z_NULL;
	if (deflateInit2(&stream_, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
	    16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
		log_fatal("Error initializing gzip encoder: %s", stream_.msg);
}

GZipEncoder::~GZipEncoder()
{
	sync();
	deflateEnd(&stream_);
}

int GZipEncoder::encode(bool flush)
{
	int ret = deflate(&stream_, flush ? Z_FINISH : Z_NO_FLUSH);
	if (ret == Z_STREAM_ERROR) {
		log_error("Error running gzip encoder: %s", stream_.msg);
		return ret;
	}
	return 0;
}
#endif

#ifdef BZIP2_FOUND
BZip2Decoder::BZip2Decoder(const std::string& path, size_t size) : Decoder(path, size)
{
	stream_.bzalloc = nullptr;
	stream_.bzfree = nullptr;
	stream_.opaque = nullptr;
	stream_.avail_in = 0;
	stream_.next_in = nullptr;
	if (BZ2_bzDecompressInit(&stream_, 0, 0) != BZ_OK)
		log_fatal("Error initializing bzip2 decoder");
}

BZip2Decoder::~BZip2Decoder()
{
	BZ2_bzDecompressEnd(&stream_);
}

int BZip2Decoder::decode()
{
	int ret = BZ2_bzDecompress(&stream_);
	if (ret != BZ_OK && ret != BZ_STREAM_END) {
		log_error("Error running bzip2 decoder");
		return ret;
	}
	return 0;
}

BZip2Encoder::BZip2Encoder(const std::string& path, size_t size) : Encoder(path, size)
{
	stream_.bzalloc = nullptr;
	stream_.bzfree = nullptr;
	stream_.opaque = nullptr;
	if (BZ2_bzCompressInit(&stream_, 9, 0, 0) != BZ_OK)
		log_fatal("Error initializing bzip2 encoder");
}

BZip2Encoder::~BZip2Encoder()
{
	sync();
	BZ2_bzCompressEnd(&stream_);
}

int BZip2Encoder::encode(bool flush)
{
	int ret = BZ2_bzCompress(&stream_, flush ? BZ_FINISH : BZ_RUN);
	if (ret == BZ_SEQUENCE_ERROR) {
		log_error("Error running bzip2 encoder");
		return ret;
	}
	return 0;
}
#endif

#ifdef LZMA_FOUND
LZMADecoder::LZMADecoder(const std::string& path, size_t size) : Decoder(path, size)
{
	stream_ = LZMA_STREAM_INIT;
	lzma_ret ret = lzma_stream_decoder(&stream_, UINT64_MAX,
	     LZMA_CONCATENATED);
	if (ret != LZMA_OK)
		log_fatal("Error initializing LZMA decoder.");
}

LZMADecoder::~LZMADecoder()
{
	lzma_end(&stream_);
}

int LZMADecoder::decode()
{
	lzma_ret ret = lzma_code(&stream_, LZMA_RUN);
	if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
		log_error("Error running LZMA decoder");
		return ret;
	}
	return 0;
}

LZMAEncoder::LZMAEncoder(const std::string& path, size_t size) : Encoder(path, size)
{
	stream_ = LZMA_STREAM_INIT;
	lzma_ret ret = lzma_easy_encoder(&stream_, LZMA_PRESET_DEFAULT,
	    LZMA_CHECK_CRC64);
	if (ret != LZMA_OK)
		log_fatal("Error initializing LZMA encoder.");
}

LZMAEncoder::~LZMAEncoder()
{
	sync();
	lzma_end(&stream_);
}

int LZMAEncoder::encode(bool flush)
{
	lzma_action action = flush ? LZMA_FINISH : LZMA_RUN;
	lzma_ret ret = lzma_code(&stream_, action);
	if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
		log_error("Error running LZMA encoder");
		return ret;
	}
	return 0;
}
#endif
