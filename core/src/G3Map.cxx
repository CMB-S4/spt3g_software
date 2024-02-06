#include <G3Map.h>

#include <pybindings.h>
#include <container_pybindings.h>
#include <serialization.h>

/* Special load/save for int64_t. */

static
int bit_count(std::map<std::string, int64_t> const &d) {
	// Returns the smallest number N such that all ints in the
	// vector could be safely expressed as intN_t.  Assumes two's
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

template <class A, typename FROM_TYPE, typename TO_TYPE>
void load_as(A &ar, std::map<std::string, TO_TYPE> &dest) {
	std::map<std::string, FROM_TYPE> temp;
	ar & cereal::make_nvp("map", temp);
	dest.insert(temp.begin(), temp.end());
}

template <class A, typename FROM_TYPE, typename TO_TYPE>
void save_as(A &ar, const std::map<std::string, FROM_TYPE> &src) {
	std::map<std::string, TO_TYPE> temp(src.begin(), src.end());
	ar & cereal::make_nvp("map", temp);
}

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
		load_as<A, int32_t, int64_t>(ar, *this);
		break;
	case 16:
		load_as<A, int16_t, int64_t>(ar, *this);
		break;
	case 8:
		load_as<A, int8_t, int64_t>(ar, *this);
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
		save_as<A, int64_t, int8_t>(ar, *this);
		break;
	case 16:
		save_as<A, int64_t, int16_t>(ar, *this);
		break;
	case 32:
		save_as<A, int64_t, int32_t>(ar, *this);
		break;
	default:
		ar & cereal::make_nvp("map",
				      cereal::base_class<std::map<std::string, int64_t> >(this));
	}
}

static
int bit_count_vector(std::map<std::string, std::vector<int64_t> > const &d) {
	// Returns the smallest number N such that all ints in the
	// vector could be safely expressed as intN_t.  Assumes two's
	// complement integers.  Return value will be between 1 and
	// 64.
	uint64_t mask = 0;
	for (auto c: d) {
		for (auto cc: c.second) {
			if (cc < 0)
				mask |= ~cc;
			else
				mask |= cc;
		}
	}
	for (int i=1; i<64; i++) {
		if (mask == 0)
			return i;
		mask >>= 1;
	}
	return 64;
}

template <class A, typename FROM_TYPE, typename TO_TYPE>
void load_vector_as(A &ar, std::map<std::string, std::vector<TO_TYPE> > &dest) {
	std::map<std::string, std::vector<FROM_TYPE> > temp;
	ar & cereal::make_nvp("map", temp);
	for (auto e: temp) {
		std::vector<TO_TYPE> v(e.second.begin(), e.second.end());
		dest[e.first] = v;
	}
}

template <class A, typename FROM_TYPE, typename TO_TYPE>
void save_vector_as(A &ar, const std::map<std::string, std::vector<FROM_TYPE> > &src) {
	std::map<std::string, std::vector<TO_TYPE> > temp;
	for (auto e: src) {
		std::vector<TO_TYPE> v(e.second.begin(), e.second.end());
		temp[e.first] = v;
	}
	ar & cereal::make_nvp("map", temp);
}

template <>
template <class A>
void G3Map<std::string, std::vector<int64_t> >::load(A &ar, const unsigned v)
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
				      cereal::base_class<std::map<std::string, std::vector<int64_t> > >(this));
		break;
	case 32:
		load_vector_as<A, int32_t, int64_t>(ar, *this);
		break;
	case 16:
		load_vector_as<A, int16_t, int64_t>(ar, *this);
		break;
	case 8:
		load_vector_as<A, int8_t, int64_t>(ar, *this);
		break;
	}
}

template <>
template <class A>
void G3Map<std::string, std::vector<int64_t> >::save(A &ar, const unsigned v) const
{
	// v == 2
	ar & cereal::make_nvp("G3FrameObject",
			      cereal::base_class<G3FrameObject>(this));
	// Count the interesting bits, and convert to nearest power of 2.
	int sig_bits = bit_count_vector(*this);
	int store_bits = 8;
	while (store_bits < sig_bits)
		store_bits *= 2;
	ar & cereal::make_nvp("store_bits", store_bits);
	switch(store_bits) {
	case 8:
		save_vector_as<A, int64_t, int8_t>(ar, *this);
		break;
	case 16:
		save_vector_as<A, int64_t, int16_t>(ar, *this);
		break;
	case 32:
		save_vector_as<A, int64_t, int32_t>(ar, *this);
		break;
	default:
		ar & cereal::make_nvp("map",
				      cereal::base_class<std::map<std::string, std::vector<int64_t> > >(this));
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
		boost::iostreams::stream<
		    boost::iostreams::back_insert_device<std::vector<char> > >
		    os(buffer);
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
	ar >> cereal::make_nvp("G3FrameObject",
	    cereal::base_class<G3FrameObject>(this));

	uint32_t len;
	ar >> cereal::make_nvp("len", len);
	for (uint32_t i = 0; i < len; i++) {
		std::pair<std::string, G3FrameObjectPtr> item;
		std::vector<char> buffer;

		ar >> cereal::make_nvp("key", item.first);
		ar >> cereal::make_nvp("value", buffer);

		boost::iostreams::array_source src((char *)&buffer[0],
		    buffer.size());
		boost::iostreams::filtering_istream fis(src);

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
G3_SERIALIZABLE_CODE(G3MapQuat);
G3_SERIALIZABLE_CODE(G3MapVectorBool);
G3_SERIALIZABLE_CODE(G3MapVectorDouble);
G3_SERIALIZABLE_CODE(G3MapVectorString);
G3_SERIALIZABLE_CODE(G3MapVectorVectorString);
G3_SERIALIZABLE_CODE(G3MapVectorComplexDouble);
G3_SERIALIZABLE_CODE(G3MapVectorTime);
G3_SERIALIZABLE_CODE(G3MapVectorQuat);

G3_SPLIT_SERIALIZABLE_CODE(G3MapInt);
G3_SPLIT_SERIALIZABLE_CODE(G3MapVectorInt);

G3_SPLIT_SERIALIZABLE_CODE(G3MapFrameObject);

PYBINDINGS("core") {
	using namespace boost::python;

	register_g3map<G3MapDouble>("G3MapDouble", "Mapping from strings to "
	    "floats");
	register_g3map<G3MapMapDouble>("G3MapMapDouble", "Mapping from strings "
	    "to maps of strings to floats. For example, "
	    "m['Det1']['Det2'] = 5.3");
	register_g3map<G3MapInt>("G3MapInt", "Mapping from strings to ints.");
	register_g3map<G3MapString>("G3MapString", "Mapping from strings to "
	    "strings.");
	register_g3map<G3MapQuat>("G3MapQuat", "Mapping from strings to "
	    "quaternions.");
	register_g3map<G3MapVectorBool>("G3MapVectorBool", "Mapping from "
	    "strings to arrays of booleans.");
	register_g3map<G3MapVectorDouble>("G3MapVectorDouble", "Mapping from "
	    "strings to arrays of floats.");
	register_g3map<G3MapVectorComplexDouble>("G3MapVectorComplexDouble",
	    "Mapping from strings to arrays of complex numbers.");
	register_g3map<G3MapVectorInt>("G3MapVectorInt", "Mapping from "
	    "strings to arrays of integers.");
	register_g3map<G3MapVectorString>("G3MapVectorString", "Mapping from "
	    "strings to lists of strings.");
	register_g3map<G3MapVectorVectorString>("G3MapVectorVectorString",
	    "Mapping from strings to lists of lists of strings.");
	register_g3map<G3MapVectorTime>("G3MapVectorTime", "Mapping from "
	    "strings to lists of G3 time objects.");
	register_g3map<G3MapVectorQuat>("G3MapVectorQuat", "Mapping from "
	    "strings to lists of quaternions.");

	// Special handling to get the object proxying right
	register_g3map<G3MapFrameObject, true>("G3MapFrameObject", "Mapping "
	    "strings to generic frame objects. Can lead to a variety of "
	    "paradoxes; please avoid general use of this class.");
}

