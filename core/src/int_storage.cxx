#include "int_storage.h"

int bit_count(const std::vector<int64_t> &d) {
	// Returns the smallest number N such that all ints in the
	// vector could be safely expressed as intN_t.  Assumes two's
	// complement integers.  Return value will be between 1 and
	// 64.
	uint64_t mask = 0;
	for (auto c: d) {
		if (c < 0)
			mask |= ~c;
		else
			mask |= c;
	}
	for (int i=1; i<64; i++) {
		if (mask == 0)
			return i;
		mask >>= 1;
	}
	return 64;
}

int bit_count(const std::map<std::string, int64_t> &d) {
	// Returns the smallest number N such that all ints in the
	// map could be safely expressed as intN_t.  Assumes two's
	// complement integers.  Return value will be between 1 and
	// 64.
	uint64_t mask = 0;
	for (auto c: d) {
		if (c.second < 0)
			mask |= ~c.second;
		else
			mask |= c.second;
	}
	for (int i=1; i<64; i++) {
		if (mask == 0)
			return i;
		mask >>= 1;
	}
	return 64;
}
