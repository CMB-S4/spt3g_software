#ifndef _G3_SKYMAPMASK_H
#define _G3_SKYMAPMASK_H

#include <G3Frame.h>
#include <G3Logging.h>

#include <vector>

/*
 * The following implements a companion object to a G3SkyMap containing
 * booleans for whether to use (true) or ignore (false) certain pixels.
 * The pixelization is stored externally in the parent map, as is the
 * (potential) translation from 2D to 1D pixel indices.
 *
 * By default, initializes to all zeroes. If use_data = True, it will
 * interpret the input sky map both as projection information and as a
 * source of data, initializing the mask to True where the input map is
 * non-zero.
 */

class G3SkyMap;

class G3SkyMapMask : public G3FrameObject {
public:
	G3SkyMapMask(const G3SkyMap &parent, bool use_data = false,
	    bool zero_nans = false, bool zero_infs = false);
	G3SkyMapMask(const G3SkyMapMask &);
	virtual ~G3SkyMapMask() {};

	// Copy
	std::shared_ptr<G3SkyMapMask> Clone(bool copy_data = true) const;

	// Return a (modifiable) pixel value
	std::vector<bool>::reference operator [] (size_t i);
	bool operator [] (size_t i) const { return this->at(i); };
	bool at(size_t i) const;

	// Logic operations:

	G3SkyMapMask &operator |=(const G3SkyMapMask &rhs);
	G3SkyMapMask &operator &=(const G3SkyMapMask &rhs);
	G3SkyMapMask &operator ^=(const G3SkyMapMask &rhs);
	G3SkyMapMask &invert(); // Basically ~=
	bool all() const;
	bool any() const;
	size_t sum() const;
	size_t size() const;
	std::vector<uint64_t> NonZeroPixels() const;

	G3SkyMapMask operator ~() const;
	G3SkyMapMask operator |(const G3SkyMapMask &rhs) const;
	G3SkyMapMask operator &(const G3SkyMapMask &rhs) const;
	G3SkyMapMask operator ^(const G3SkyMapMask &rhs) const;
	G3SkyMapMask operator ==(const G3SkyMapMask &rhs) const;
	G3SkyMapMask operator !=(const G3SkyMapMask &rhs) const;

	void ApplyMask(const G3SkyMapMask &rhs, bool inverse=false);

	// Information
	bool IsCompatible(const G3SkyMap &map) const;
	bool IsCompatible(const G3SkyMapMask &mask) const;
	bool IsDense() const { return true; }

	// The map for projection info
	std::shared_ptr<const G3SkyMap> Parent() const { return parent_; }

	// Return a 1-or-0 sky-map with the contents of the mask (e.g. for plotting)
	std::shared_ptr<G3SkyMap> MakeBinaryMap() const;

	class const_iterator {
	public:
		typedef std::pair<uint64_t, bool> value_type;
		typedef value_type & reference;
		typedef value_type * pointer;

		const_iterator(const G3SkyMapMask &mask, bool begin);

		bool operator==(const const_iterator & other) const {
			return index_ == other.index_;
		};

		bool operator!=(const const_iterator & other) const {
			return index_ != other.index_;
		};

		reference operator*() { return value_; };
		pointer operator->() { return &value_; };

		const_iterator operator++();
		const_iterator operator++(int) { const_iterator i = *this; ++(*this); return i; }

	private:
		uint64_t index_;
		value_type value_;
		const G3SkyMapMask &mask_;

		void set_value() {
			value_.first = index_;
			value_.second = mask_.at(index_);
		}
	};

	const_iterator begin() const { return const_iterator(*this, true); };
	const_iterator end() const { return const_iterator(*this, false); };

private:
	G3SkyMapMask() {} // Fake out for serialization
	template <class A> void load(A &ar, const unsigned v);
	template <class A> void save(A &ar, const unsigned v) const;
	friend class cereal::access;

	std::vector<bool> data_;
	std::shared_ptr<const G3SkyMap> parent_;

	SET_LOGGER("G3SkyMapMask");
};

namespace cereal {
  template <class A> struct specialize<A, G3SkyMapMask, cereal::specialization::member_load_save> {};
}

G3_POINTERS(G3SkyMapMask);
G3_SERIALIZABLE(G3SkyMapMask, 2);

#endif

