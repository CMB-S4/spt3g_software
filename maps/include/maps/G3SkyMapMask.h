#ifndef _G3_SKYMAP_H
#define _G3_SKYMAP_H

#include <G3SkyMap.h>
#include <vector>

#include <boost/python.hpp>

/*
 * The following implements a companion object to a G3SkyMap containing
 * booleans for whether to use (true) or ignore (false) certain pixels.
 * The pixelization is stored externally in the parent map, as is the
 * (potential) translation from 2D to 1D pixel indices.
 */

class G3SkyMapMask : public G3FrameObject {
public:
	G3SkyMapMask(const G3SkyMap &parent);
	virtual ~G3SkyMapMask() {};

	// Return a (modifiable) pixel value
	bool &operator [] (size_t i);
	bool operator [] (size_t i) const { return this->at(i); };
	bool at(size_t i) const;

	// Logic operations:

	G3SkyMapMask &operator |=(const G3SkyMapMask &rhs);
	G3SkyMapMask &operator &=(const G3SkyMapMask &rhs);
	G3SkyMapMask &operator ^=(const G3SkyMapMask &rhs);
	G3SkyMapMask &invert(); // Basically ~=

	G3SkyMapMask operator ~();
	G3SkyMapMask operator |(const G3SkyMapMask &rhs);
	G3SkyMapMask operator &(const G3SkyMapMask &rhs);
	G3SkyMapMask operator ^(const G3SkyMapMask &rhs);

	bool IsDense() const { return true; }

	// The map for projection info
	G3SkyMapConstPtr Parent() const;

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

