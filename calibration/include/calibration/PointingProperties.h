#ifndef POINTINGPROPERTIES_H
#define POINTINGPROPERTIES_H

#include <G3Frame.h>
#include <G3Map.h>
#include <string>

/*
 * Class containing offline pointing model parameters.
 */
class PointingProperties : public G3FrameObject
{
public:
	PointingProperties() : tiltLat(NAN), tiltHA(NAN), tiltMag(NAN),
	    tiltAngle(NAN) {}
	std::string Description() const;

	double tiltLat;
	double tiltHA;
	double tiltMag;
	double tiltAngle;

	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(PointingProperties);
G3_SERIALIZABLE(PointingProperties, 4);

G3MAP_OF(std::string, PointingProperties, PointingPropertiesMap);

#endif

