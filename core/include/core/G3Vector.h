#ifndef _G3_VECTOR_H
#define _G3_VECTOR_H

#include <G3Frame.h>
#include <G3TimeStamp.h>
#include <vector>
#include <sstream>
#include <complex>

#include <cereal/types/complex.hpp>
#include <cereal/types/vector.hpp>

/*
 * Generic container that is an std::vector of things and is also a frame object
 */

template <typename Value>
class G3Vector : public G3FrameObject, public std::vector<Value> {
public:
	G3Vector() {}
	G3Vector(typename std::vector<Value>::size_type s) :
	    std::vector<Value>(s) {}
	G3Vector(typename std::vector<Value>::size_type s, const Value &val) :
	    std::vector<Value>(s, val) {}
	G3Vector(const G3Vector &r) : std::vector<Value>(r) {}
	template <typename Iterator> G3Vector(Iterator l, Iterator r) :
	    std::vector<Value>(l, r) {}

	template <class A> void serialize(A &ar, const unsigned v)
	{
		ar & cereal::make_nvp("G3FrameObject",
		    cereal::base_class<G3FrameObject>(this));
		ar & cereal::make_nvp("vector",
		    cereal::base_class<std::vector<Value> >(this));
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
		s << "[";
		if (this->size() == 1) {
			s << (*this)[0];
		} else if (this->size() > 1) {
			for (int i = 0; i < this->size()-1; i++)
				s << (*this)[i] << ", ";
			s << (*this)[this->size() - 1];
		}
		s << "]";
		return s.str();
	}
};

#define G3VECTOR_OF(x, y) \
typedef G3Vector< x > y; \
namespace cereal { \
	template <class A> struct specialize<A, y, cereal::specialization::member_serialize> {}; \
} \
G3_POINTERS(y); \
G3_SERIALIZABLE(y, 1);

G3VECTOR_OF(std::complex<double>, G3VectorComplexDouble);
G3VECTOR_OF(double, G3VectorDouble);
G3VECTOR_OF(int32_t, G3VectorInt);
G3VECTOR_OF(uint8_t, G3VectorUnsignedChar);
G3VECTOR_OF(std::string, G3VectorString);
G3VECTOR_OF(G3VectorString, G3VectorVectorString);
G3VECTOR_OF(G3FrameObjectPtr, G3VectorFrameObject);
G3VECTOR_OF(G3Time, G3VectorTime);

#endif

