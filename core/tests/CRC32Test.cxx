#include <G3Test.h>

#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

extern "C" unsigned long crc32c(unsigned long crc, const uint8_t *buf, unsigned int len);
extern "C" const char *crc32c_impl(void);

TEST_GROUP(CRC32C)

// Bitwise reference implementation, with no pre- or post-inversion to match
// the convention used in G3 files
static uint32_t crc32c_ref(uint32_t crc, const uint8_t *buf, size_t len)
{
	for (size_t i = 0; i < len; i++) {
		crc ^= buf[i];
		for (int k = 0; k < 8; k++)
			crc = (crc & 1) ? (crc >> 1) ^ 0x82f63b78 : crc >> 1;
	}
	return crc;
}

static std::vector<uint8_t> test_data(size_t len)
{
	std::vector<uint8_t> data(len);
	uint32_t x = 12345;
	for (size_t i = 0; i < len; i++) {
		x = x * 1103515245 + 12345;
		data[i] = x >> 24;
	}
	return data;
}

TEST(KnownValue){
	const uint8_t check[] = "123456789";
	ENSURE_EQUAL(crc32c(0, check, 9), 0x58e3fa20ul);
	ENSURE_EQUAL(crc32c(0, check, 0), 0ul);
}

TEST(LengthsAndAlignments){
	// Lengths around the boundaries of the parallel-stream blocks used by
	// the hardware implementation (3 x 256 and 3 x 8192 bytes)
	const size_t lengths[] = {1, 7, 8, 9, 63, 767, 768, 769, 1000,
	    24575, 24576, 24577, 3 * 24576 + 3 * 256 + 13, 100003};
	std::vector<uint8_t> data = test_data(100003 + 8);

	for (size_t len : lengths) {
		for (size_t offset = 0; offset < 8; offset++) {
			const uint8_t *buf = data.data() + offset;
			ENSURE_EQUAL(crc32c(0, buf, len),
			    (unsigned long)crc32c_ref(0, buf, len),
			    "length " + std::to_string(len) + ", offset " +
			    std::to_string(offset));
		}
	}
}

TEST(Chained){
	// CRCs of consecutive pieces of a buffer must match the CRC of the whole
	std::vector<uint8_t> data = test_data(100000);
	uint32_t expected = crc32c_ref(0, data.data(), data.size());

	const size_t splits[] = {1, 5, 768, 24576, 50001, 99999};
	for (size_t split : splits) {
		unsigned long crc = crc32c(0, data.data(), split);
		crc = crc32c(crc, data.data() + split, data.size() - split);
		ENSURE_EQUAL(crc, (unsigned long)expected,
		    "split at " + std::to_string(split));
	}
}

// Throughput in GiB/s of the best of several passes over a buffer
static double crc32c_throughput(const std::vector<uint8_t> &data, int passes)
{
	double best = 1e30;
	volatile unsigned long sink = 0;
	for (int i = 0; i < passes; i++) {
		auto start = std::chrono::steady_clock::now();
		sink = crc32c(sink, data.data(), data.size());
		std::chrono::duration<double> dt =
		    std::chrono::steady_clock::now() - start;
		if (dt.count() < best)
			best = dt.count();
	}
	return data.size() / best / (1 << 30);
}

// Report which implementation is in use and how fast it is. Output is only
// shown when running in verbose mode (core-test --verbose).
TEST(Diagnostics){
#if defined(__x86_64__)
	const char *arch = "x86_64";
#elif defined(__aarch64__)
	const char *arch = "aarch64";
#elif defined(__i386__)
	const char *arch = "i386";
#else
	const char *arch = "other";
#endif
	double cached = crc32c_throughput(test_data(256 << 10), 200);
	double large = crc32c_throughput(test_data(64 << 20), 3);

	std::cout << std::fixed << std::setprecision(1) << "\n"
	    << "  Target architecture: " << arch << "\n"
	    << "  CRC32C implementation: " << crc32c_impl() << "\n"
	    << "  Throughput, 256 KiB buffer (in cache): " << cached << " GiB/s\n"
	    << "  Throughput, 64 MiB buffer: " << large << " GiB/s\n";
}
