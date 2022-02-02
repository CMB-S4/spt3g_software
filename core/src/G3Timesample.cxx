#include <pybindings.h>

#include <numeric>
#include <algorithm>
#include <iostream>
#include <boost/python.hpp>

#include <container_pybindings.h>
#include <G3Timesample.h>


// Templates for vector operations.  These work on any std::vector,
// and thus on any G3Vector.

// Returns the vector size if the type matches, and -1 otherwise.
template <typename g3vectype>
inline
ssize_t vect_size(const G3FrameObjectPtr &vp)
{
	auto v = boost::dynamic_pointer_cast<const g3vectype>(vp);
	return v == nullptr ? -1 : v->size();
}

// Concatenates two compatible vectors.
template <typename vectype>
inline
void vect_concat(vectype &dest, const vectype &src1, const vectype &src2)
{
	dest.clear();
	dest.reserve(src1.size() + src2.size());
	dest.insert(dest.end(), src1.begin(), src1.end());
	dest.insert(dest.end(), src2.begin(), src2.end());
}

// This returns the length of the vector if it's a valid type (within
// the context of G3Timesample).  For invalid types, it returns -1.
inline
ssize_t g3_vect_test_and_size(const G3FrameObjectPtr &vp)
{
	ssize_t vsize = -1;
	if ((vsize = vect_size<G3VectorDouble>(vp)) < 0 &&
	    (vsize = vect_size<G3VectorInt   >(vp)) < 0 &&
	    (vsize = vect_size<G3VectorBool  >(vp)) < 0 &&
	    (vsize = vect_size<G3VectorString>(vp)) < 0)
		return -1;
	return vsize;
}

// Concatenates two compatible vectors and returns the result.
// Returns nullptr on any incompatibility.
template <typename g3vectype>
inline
G3FrameObjectPtr test_and_concat(const G3FrameObjectPtr &src1,
				 const G3FrameObjectPtr &src2)
{
	auto v1 = boost::dynamic_pointer_cast<const g3vectype>(src1);
	auto v2 = boost::dynamic_pointer_cast<const g3vectype>(src2);
	if (v1 == nullptr || v2 == nullptr)
		return nullptr;
	auto output = boost::shared_ptr<g3vectype>(new g3vectype());
	vect_concat(*output, *v1, *v2);
	return output;
}

// Reorders elements of a vector.
template <typename vectype>
inline
void vect_reorder(vectype &dest, vectype &src, const std::vector<size_t> &idx)
{
	dest.clear();
	dest.resize(src.size());
	assert(src.size() == idx.size());
	for (size_t i = 0; i < idx.size(); i++)
		dest[i] = src[idx[i]];
}

// Reorders elements of a vector, in place.
// Returns false on any incompatibility.
template <typename g3vectype>
inline
bool test_and_reorder(G3FrameObjectPtr &src, const std::vector<size_t> &idx)
{
	auto v = boost::dynamic_pointer_cast<g3vectype>(src);
	if (v == nullptr)
		return false;
	g3vectype v2 = *v; // copy
	vect_reorder(*v, v2, idx);
	return true;
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

bool G3TimesampleMap::Check() const
{
	int n = times.size();

	for (auto item = begin(); item != end(); ++item) {
		auto name = item->first;
		auto el = item->second;

		int check_len;
		if ((check_len = g3_vect_test_and_size(el)) < 0) {
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

G3TimesampleMap G3TimesampleMap::Concatenate(const G3TimesampleMap &other) const
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
	vect_concat(output.times, times, other.times);

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
			(catted = test_and_concat<G3VectorBool>(
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

void G3TimesampleMap::Sort()
{
	Check();

	// Short circuit no-op
	if (std::is_sorted(times.begin(), times.end()))
		return;

	// Get time vector sort order
	std::vector<size_t> idx(times.size());
	std::iota(idx.begin(), idx.end(), 0);
	std::stable_sort(idx.begin(), idx.end(),
	    [this](size_t i1, size_t i2) { return times[i1] < times[i2]; });

	// Sort all the things
	G3VectorTime oldtimes = times;
	vect_reorder(times, oldtimes, idx);

	for (auto item = begin(); item != end(); item++) {
		if (
			!test_and_reorder<G3VectorDouble>(item->second, idx) &&
			!test_and_reorder<G3VectorInt>(item->second, idx) &&
			!test_and_reorder<G3VectorBool>(item->second, idx) &&
			!test_and_reorder<G3VectorString>(item->second, idx)
			) {
			std::ostringstream s;
			s << "Vector type not supported for key: " << item->first << "\n";
			throw g3timesample_exception(s.str());
		}
	}
}


// Safety-ized for python.

static
void safe_set_item(G3TimesampleMap &self, const std::string key,
		   G3FrameObjectPtr value)
{
	int check_len = g3_vect_test_and_size(value);
	if (check_len < 0) {
		std::ostringstream s;
		s << "Cannot add member (" << key << "): "
		  << "not a supported vector type.";
		throw g3timesample_exception(s.str());
	}
	if (check_len != self.times.size()) {
		std::ostringstream s;
		s << "Cannot add member (" << key << "): "
		  << "not the same length as .times.";
		throw g3timesample_exception(s.str());
	}
	self[key] = value;
}

static
void safe_set_times(G3TimesampleMap &self, G3VectorTime _times)
{
	// Only allow this if it doesn't upset consistency.  We will
	// assume that, coming in, we're internally consistent.
	if (_times.size() != self.times.size() && self.size() != 0) {
		std::ostringstream s;
		s << "Cannot set .times because it conflicts with "
		  << "the established number of samples (" << self.times.size()
		  << ").";
		throw g3timesample_exception(s.str());
	}
	self.times = _times;
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
	.def("__setitem__", &safe_set_item)
	.def_pickle(g3frameobject_picklesuite<G3TimesampleMap>())
	// Extensions for G3TimesampleMap are here:
	.add_property("times", &G3TimesampleMap::times, &safe_set_times,
	  "Times vector.  Setting this stores a copy, but getting returns a reference.")
	.def("check", &G3TimesampleMap::Check, "Check for internal "
          "consistency.  Raises ValueError if there are problems.")
	.def("concatenate", &G3TimesampleMap::Concatenate,
          "Concatenate two compatible G3TimesampleMap.")
	.def("sort", &G3TimesampleMap::Sort,
          "Sort all element vectors by time, in-place.")
	;
	register_pointer_conversions<G3TimesampleMap>();

	bp::register_exception_translator<g3timesample_exception>(&translate_ValueError);
}
