#include <G3Map.h>

#include <pybindings.h>
#include <container_pybindings.h>
#include <serialization.h>

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

G3_SERIALIZABLE_CODE(G3MapInt);
G3_SERIALIZABLE_CODE(G3MapDouble);
G3_SERIALIZABLE_CODE(G3MapMapDouble);
G3_SERIALIZABLE_CODE(G3MapString);
G3_SERIALIZABLE_CODE(G3MapVectorInt);
G3_SERIALIZABLE_CODE(G3MapVectorDouble);
G3_SERIALIZABLE_CODE(G3MapVectorString);
G3_SERIALIZABLE_CODE(G3MapVectorVectorString);
G3_SERIALIZABLE_CODE(G3MapVectorComplexDouble);
G3_SERIALIZABLE_CODE(G3MapVectorTime);

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

	// Special handling to get the object proxying right
	register_g3map<G3MapFrameObject, true>("G3MapFrameObject", "Mapping "
	    "strings to generic frame objects. Can lead to a variety of "
	    "paradoxes; please avoid general use of this class.");
}

