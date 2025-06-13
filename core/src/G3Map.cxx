#include <G3Map.h>

#include <pybindings.h>
#include <container_pybindings.h>
#include <serialization.h>
#include "int_storage.h"

/* Special load/save for int64_t, using the same encoding at G3VectorInt */

template <>
template <class A>
void G3Map<std::string, int64_t>::load(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	int store_bits = 32;
	if (v >= 2)
		ar & cereal::make_nvp("store_bits", store_bits);

	switch(store_bits) {
	case 64:
		ar & cereal::make_nvp("map",
		    cereal::base_class<std::map<std::string, int64_t> >(this));
		break;
	case 32:
		load_as<A, int32_t>(ar, *this);
		break;
	case 16:
		load_as<A, int16_t>(ar, *this);
		break;
	case 8:
		load_as<A, int8_t>(ar, *this);
		break;
	}
}

template <>
template <class A>
void G3Map<std::string, int64_t>::save(A &ar, const unsigned v) const
{
	// v == 2
	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));
	// Count the interesting bits, and convert to nearest power of 2.
	int sig_bits = bit_count(*this);
	int store_bits = 8;
	while (store_bits < sig_bits)
		store_bits *= 2;
	ar & cereal::make_nvp("store_bits", store_bits);
	switch(store_bits) {
	case 8:
		save_as<A, int8_t>(ar, *this);
		break;
	case 16:
		save_as<A, int16_t>(ar, *this);
		break;
	case 32:
		save_as<A, int32_t>(ar, *this);
		break;
	default:
		ar & cereal::make_nvp("map",
		    cereal::base_class<std::map<std::string, int64_t> >(this));
	}
}

template <>
template <class A>
void G3Map<std::string, std::vector<int64_t> >::load(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));

	if (v == 1) {
		std::map<std::string, std::vector<int32_t> > temp;
		ar & cereal::make_nvp("map", temp);
		for (auto const &i: temp) {
			std::vector<int64_t> v(i.second.begin(), i.second.end());
			(*this)[i.first] = v;
		}

		return;
	}

	uint32_t len;
	ar & cereal::make_nvp("len", len);

	for (uint32_t i = 0; i < len; i++) {
		std::pair<std::string, std::vector<int64_t> > item;
		ar & cereal::make_nvp("key", item.first);
		int store_bits;
		ar & cereal::make_nvp("store_bits", store_bits);

		switch(store_bits) {
		case 64:
			ar & cereal::make_nvp("vector", item.second);
			break;
		case 32:
			load_as<A, int32_t>(ar, item.second);
			break;
		case 16:
			load_as<A, int16_t>(ar, item.second);
			break;
		case 8:
			load_as<A, int8_t>(ar, item.second);
			break;
		}

		this->insert(item);
	}
}

template <>
template <class A>
void G3Map<std::string, std::vector<int64_t> >::save(A &ar, const unsigned v) const
{
	// v == 2
	ar & cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));

	uint32_t len = size();
	ar & cereal::make_nvp("len", len);

	for (auto const &i: *this) {
		ar & cereal::make_nvp("key", i.first);

		// Count the interesting bits, and convert to nearest power of 2.
		int sig_bits = bit_count(i.second);
		int store_bits = 8;
		while (store_bits < sig_bits)
			store_bits *= 2;
		ar & cereal::make_nvp("store_bits", store_bits);

		switch(store_bits) {
		case 8:
			save_as<A, int8_t>(ar, i.second);
			break;
		case 16:
			save_as<A, int16_t>(ar, i.second);
			break;
		case 32:
			save_as<A, int32_t>(ar, i.second);
			break;
		default:
			ar & cereal::make_nvp("vector", i.second);
		}
	}
}

template <class A> void G3MapFrameObject::save(A &ar, const unsigned v) const
{
	ar << cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));

	uint32_t len = size();
	ar << cereal::make_nvp("len", len);

	for (auto i = begin(); i != end(); i++) {
		ar << cereal::make_nvp("key", i->first);

		std::vector<char> buffer;
		G3BufferOutputStream os(buffer);
		{
			A subar(os);
			subar << cereal::make_nvp("item",
			    i->second);
		}
		os.flush();
		ar << cereal::make_nvp("value", buffer);
	}
}

template <class A> void G3MapFrameObject::load(A &ar, const unsigned v)
{
	G3_CHECK_VERSION(v);

	ar >> cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));

	uint32_t len;
	ar >> cereal::make_nvp("len", len);
	for (uint32_t i = 0; i < len; i++) {
		std::pair<std::string, G3FrameObjectPtr> item;
		std::vector<char> buffer;

		ar >> cereal::make_nvp("key", item.first);
		ar >> cereal::make_nvp("value", buffer);

		G3BufferInputStream fis(buffer);
		A subar(fis);
		subar >> cereal::make_nvp("item", item.second);
		this->insert(item);
	}
}

std::string G3MapFrameObject::Summary() const
{
	std::ostringstream s;
	s << this->size() << " elements";
	return s.str();
}

std::string G3MapFrameObject::Description() const
{
	std::ostringstream s;
	s << '{';
	for (auto i = this->begin(); i != this->end(); i++)
		s << i->first << ": " << i->second->Summary() << ", ";
	s << '}';
	return s.str();
}

G3_SERIALIZABLE_CODE(G3MapDouble);
G3_SERIALIZABLE_CODE(G3MapMapDouble);
G3_SERIALIZABLE_CODE(G3MapString);
G3_SERIALIZABLE_CODE(G3MapTime);
G3_SERIALIZABLE_CODE(G3MapVectorBool);
G3_SERIALIZABLE_CODE(G3MapVectorDouble);
G3_SERIALIZABLE_CODE(G3MapVectorString);
G3_SERIALIZABLE_CODE(G3MapVectorVectorString);
G3_SERIALIZABLE_CODE(G3MapVectorComplexDouble);
G3_SERIALIZABLE_CODE(G3MapVectorTime);

G3_SPLIT_SERIALIZABLE_CODE(G3MapInt);
G3_SPLIT_SERIALIZABLE_CODE(G3MapVectorInt);

G3_SPLIT_SERIALIZABLE_CODE(G3MapFrameObject);

PYBINDINGS("core", scope) {
	register_g3map<G3MapDouble>(scope, "G3MapDouble", "Mapping from strings to "
	    "floats");
	register_g3map<G3MapMapDouble>(scope, "G3MapMapDouble", "Mapping from strings "
	    "to maps of strings to floats. For example, "
	    "m['Det1']['Det2'] = 5.3");
	register_g3map<G3MapInt>(scope, "G3MapInt", "Mapping from strings to ints.");
	register_g3map<G3MapString>(scope, "G3MapString", "Mapping from strings to "
	    "strings.");
	register_g3map<G3MapTime>(scope, "G3MapTime", "Mapping from strings to "
	    "G3Times.");
	register_g3map<G3MapVectorBool>(scope, "G3MapVectorBool", "Mapping from "
	    "strings to arrays of booleans.");
	register_g3map<G3MapVectorDouble>(scope, "G3MapVectorDouble", "Mapping from "
	    "strings to arrays of floats.");
	register_g3map<G3MapVectorComplexDouble>(scope, "G3MapVectorComplexDouble",
	    "Mapping from strings to arrays of complex numbers.");
	register_g3map<G3MapVectorInt>(scope, "G3MapVectorInt", "Mapping from "
	    "strings to arrays of integers.");
	register_g3map<G3MapVectorString>(scope, "G3MapVectorString", "Mapping from "
	    "strings to lists of strings.");
	register_g3map<G3MapVectorVectorString>(scope, "G3MapVectorVectorString",
	    "Mapping from strings to lists of lists of strings.");
	register_g3map<G3MapVectorTime>(scope, "G3MapVectorTime", "Mapping from "
	    "strings to lists of G3 time objects.");

	register_g3map<G3MapFrameObject>(scope, "G3MapFrameObject", "Mapping "
	    "strings to generic frame objects. Can lead to a variety of "
	    "paradoxes; please avoid general use of this class.");
}

