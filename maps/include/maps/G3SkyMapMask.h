#ifndef _G3_SKYMAPMASK_H
#define _G3_SKYMAPMASK_H

#include <maps/G3SkyMap.h>
#include <vector>

#include <boost/python.hpp>

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

class G3SkyMapMask : public G3FrameObject {
public:
	G3SkyMapMask(const G3SkyMap &parent, bool use_data = false);
	G3SkyMapMask(const G3SkyMapMask &);
	virtual ~G3SkyMapMask() {};

	// Return a (modifiable) pixel value
	std::vector<bool>::reference operator [] (size_t i);
	bool operator [] (size_t i) const { return this->at(i); };
	bool at(size_t i) const;

	// Logic operations:

	G3SkyMapMask &operator |=(const G3SkyMapMask &rhs);
	G3SkyMapMask &operator &=(const G3SkyMapMask &rhs);
	G3SkyMapMask &operator ^=(const G3SkyMapMask &rhs);
	G3SkyMapMask &invert(); // Basically ~=

	G3SkyMapMask operator ~() const;
	G3SkyMapMask operator |(const G3SkyMapMask &rhs) const;
	G3SkyMapMask operator &(const G3SkyMapMask &rhs) const;
	G3SkyMapMask operator ^(const G3SkyMapMask &rhs) const;

	// Information
	bool IsCompatible(const G3SkyMap &map) const { return map.IsCompatible(*Parent()); }
	bool IsCompatible(const G3SkyMapMask &mask) const { return mask.Parent()->IsCompatible(*Parent()); }
	bool IsDense() const { return true; }

	// The map for projection info
	G3SkyMapConstPtr Parent() const { return parent_; }

	// Return a 1-or-0 sky-map with the contents of the mask (e.g. for plotting)
	G3SkyMapPtr MakeBinaryMap() const;

private:
	G3SkyMapMask() {} // Fake out for serialization
	template <class A> void serialize(A &ar, const unsigned v);
	friend class cereal::access;

	std::vector<bool> data_;
	G3SkyMapConstPtr parent_;

	SET_LOGGER("G3SkyMapMask");
};

G3_POINTERS(G3SkyMapMask);
G3_SERIALIZABLE(G3SkyMapMask, 1);

#endif

