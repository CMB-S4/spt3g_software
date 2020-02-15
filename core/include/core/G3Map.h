#ifndef _G3_MAP_H
#define _G3_MAP_H

#include <G3Frame.h>
#include <G3Vector.h>
#include <map>
#include <sstream>
#include <complex>

#include <cereal/types/complex.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>

/*
 * Generic container that is an std::map of things and is also a frame object
 */

template <typename Key, typename Value>
class G3Map : public G3FrameObject, public std::map<Key, Value> {
public:
	template <class A> void serialize(A &ar, unsigned v)
	{
		ar & cereal::make_nvp("G3FrameObject",
		    cereal::base_class<G3FrameObject>(this));
		ar & cereal::make_nvp("map",
		    cereal::base_class<std::map<Key, Value> >(this));
	}

	std::string Summary() const
	{
		if (this->size() < 5)
			return Description();

		std::ostringstream s;
		s << this->size() << " elements";
		return s.str();
	}

	std::string Description() const
	{
		std::ostringstream s;
		s << '{';
		for (auto i = this->begin(); i != this->end(); i++)
			s << i->first << ", ";
		s << '}';
		return s.str();
	}
};

// Work around a bug in Boost related to serialization of multiple types
// through the same base pointer in the same archive by providing a separate
// class
class G3MapFrameObject : public G3FrameObject, public std::map<std::string, G3FrameObjectPtr> {
public:
	template <class A> void save(A &ar, const unsigned v) const;
	template <class A> void load(A &ar, const unsigned v);

	std::string Summary() const;
	std::string Description() const;
};

#define G3MAP_OF(key, value, name) \
typedef G3Map< key, value > name; \
namespace cereal { \
        template <class A> struct specialize<A, name, cereal::specialization::member_serialize> {}; \
} \
G3_POINTERS(name); \
G3_SERIALIZABLE(name, 1);

G3MAP_OF(std::string, double, G3MapDouble);
G3MAP_OF(std::string, G3MapDouble, G3MapMapDouble);
G3MAP_OF(std::string, std::vector<double>, G3MapVectorDouble);
G3MAP_OF(std::string, std::vector<int32_t>, G3MapVectorInt);
G3MAP_OF(std::string, std::vector<std::string>, G3MapVectorString);
G3MAP_OF(std::string, G3VectorVectorString, G3MapVectorVectorString);
G3MAP_OF(std::string, std::vector<std::complex<double> >, G3MapVectorComplexDouble);
G3MAP_OF(std::string, G3VectorTime, G3MapVectorTime);
G3MAP_OF(std::string, int32_t, G3MapInt);
G3MAP_OF(std::string, std::string, G3MapString);

namespace cereal {
        template <class A> struct specialize<A, G3MapFrameObject, cereal::specialization::member_load_save> {};
}
G3_POINTERS(G3MapFrameObject);
G3_SERIALIZABLE(G3MapFrameObject, 1);
#endif

