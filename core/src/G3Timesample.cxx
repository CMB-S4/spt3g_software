#include <pybindings.h>

#include <iostream>
#include <boost/python.hpp>

#include <container_pybindings.h>
#include <G3Timesample.h>
//#include <exceptions.h>

// Templates for vector operations.  g3vectype might be, for example,
// G3VectorDouble.  One could go all the way and SFINAE this... but
// perhaps not everything should require a black belt.

template <typename g3vectype>
inline
int g3_vect_size(const G3FrameObjectPtr vp)
{
	auto v = boost::dynamic_pointer_cast<const g3vectype>(vp);
	return v == nullptr ? -1 : v->size();
}

template <typename g3vectype>
inline
void g3_concat(g3vectype &output, const g3vectype &src1, const g3vectype &src2)
{
	auto dest = output.begin();
	for (auto p1: src1)
		*(dest++) = p1;
	for (auto p2: src2)
		*(dest++) = p2;
}

template <typename g3vectype>
inline
G3FrameObjectPtr test_and_concat(const G3FrameObjectPtr src1, const G3FrameObjectPtr src2)
{
	auto v1 = boost::dynamic_pointer_cast<const g3vectype>(src1);
	auto v2 = boost::dynamic_pointer_cast<const g3vectype>(src2);
	if (v1 == nullptr || v2 == nullptr)
		return nullptr;
	boost::shared_ptr<g3vectype> outputp(new g3vectype());
	outputp->resize(v1->size() + v2->size());
	g3_concat(*outputp, *v1, *v2);
	return outputp;
}


/* G3TimesampleMap */

std::string G3TimesampleMap::Description() const
{
	std::ostringstream s;
	s << "<co-sampled vectors with " << times.size() << " samples>{";
	for (auto i = this->begin(); i != this->end(); ) {
		s << i->first;
		if (++i != this->end())
			s << ", ";
	}
	s << "}";
	return s.str();
}

std::string G3TimesampleMap::Summary() const
{
	return Description();
}

template <class A> void G3TimesampleMap::serialize(A &ar, unsigned v)
{
	using namespace cereal;
	ar & make_nvp("parent", base_class<G3MapFrameObject>(this));
	ar & make_nvp("times", times);
}

bool G3TimesampleMap::Check()
{
	int n = times.size();

	// How to polymorph?  This is how.
	for (auto item = begin(); item != end(); ++item) {
		auto name = item->first;
		auto el = item->second;

		int check_type = 0;
		int check_len = -1;

		// Try to get a length...
		if ((check_len = g3_vect_size<G3VectorDouble>(el)) < 0 &&
		    (check_len = g3_vect_size<G3VectorInt   >(el)) < 0 &&
		    (check_len = g3_vect_size<G3VectorString>(el)) < 0) {
			std::ostringstream s;
			s << "Vector type not supported for key: " << name << "\n";
			throw g3timesample_exception(s.str());
		}

		if (check_len != n) {
			std::ostringstream s;
			s << "Vector not same length as .times: " << name << "\n";
			throw g3timesample_exception(s.str());
		}
	}
	return true;
}

const
G3TimesampleMap G3TimesampleMap::Concatenate(const G3TimesampleMap other)
{
	// Check that all keys in other are in this.
	for (auto item = other.begin(); item != other.end(); ++item) {
		if (find(item->first) == end()) {
			std::ostringstream s;
			s << "Inconsistent keys; " << item->first << " on right only.";
			throw g3timesample_exception(s.str());
		}
	}

	int n_cat = times.size() + other.times.size();
	G3TimesampleMap output;
	output.times.resize(n_cat);
	g3_concat(output.times, times, other.times);

	for (auto item = begin(); item != end(); ++item) {
		auto oitem = other.find(item->first);
		if (oitem == other.end()) {
			std::ostringstream s;
			s << "Inconsistent keys; " << item->first << " on left only.";
			throw g3timesample_exception(s.str());
		}

		G3FrameObjectPtr catted;
		if (
			(catted = test_and_concat<G3VectorDouble>(
				item->second, oitem->second)) != nullptr ||
			(catted = test_and_concat<G3VectorInt>(
				item->second, oitem->second)) != nullptr ||
			(catted = test_and_concat<G3VectorString>(
				item->second, oitem->second)) != nullptr
			) {
			output.insert(std::make_pair(item->first, catted));
		} else {
			std::ostringstream s;
			s << "Vector type not supported for key: " << item->first << "\n";
			throw g3timesample_exception(s.str());
		}
	}

	return output;
}


G3_SERIALIZABLE_CODE(G3TimesampleMap);

static void translate_ValueError(g3timesample_exception const& e)
{
    PyErr_SetString(PyExc_ValueError, e.msg_for_python().c_str());
}


PYBINDINGS("core")
{
	// This is based on register_g3map macro.
	bp::class_<G3TimesampleMap, bp::bases<G3FrameObject,
	    std::map<typename G3TimesampleMap::key_type,
                     typename G3TimesampleMap::mapped_type> >,
	    boost::shared_ptr<G3TimesampleMap> >("G3TimesampleMap",
              "Mapping from string to vectors of data, with an associated "
              "vector of timestamps.  This object is for storing multiple "
              "co-sampled vectors with a single set of (irregular) timestamps.")
	.def(bp::init<const G3TimesampleMap &>())
	.def(bp::std_map_indexing_suite<G3TimesampleMap, true>())
	.def_pickle(g3frameobject_picklesuite<G3TimesampleMap>())
	// Extensions for G3TimesampleMap are here:
	.def_readwrite("times", &G3TimesampleMap::times, "Timestamp vector.")
	.def("Check", &G3TimesampleMap::Check, "Check for internal "
          "consistency.  Raises ValueError if there are problems.")
	.def("Concatenate", &G3TimesampleMap::Concatenate,
          "Concatenate two compatible G3TimesampleMap.")
	;
	register_pointer_conversions<G3TimesampleMap>();

	bp::register_exception_translator<g3timesample_exception>(&translate_ValueError);
}
