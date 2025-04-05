#ifndef _G3_SUPERTIMESTREAM_H
#define _G3_SUPERTIMESTREAM_H

#include <G3Logging.h>
#include <G3Timestream.h>

class G3SuperTimestream : public G3TimestreamMap {
private:
	G3SuperTimestream() {};

	friend class cereal::access;

	template <class A> void load(A &ar, unsigned v);
	template <class A> void save(A &ar, unsigned v) const;

	SET_LOGGER("G3SuperTimestream");
};

namespace cereal {
	template <class A> struct specialize<A, G3SuperTimestream, cereal::specialization::member_load_save> {};
}

G3_POINTERS(G3SuperTimestream);
G3_SERIALIZABLE(G3SuperTimestream, 0);

#endif
