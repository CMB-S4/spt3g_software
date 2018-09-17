#ifndef BOLOPROPERTIES_H
#define BOLOPROPERTIES_H

#include <G3Frame.h>
#include <G3Map.h>
#include <string>

/*
 * Class containing static physical properties of a bolometer (i.e. focal plane
 * position, but not gain).
 */
class BolometerProperties : public G3FrameObject
{
public:
	BolometerProperties() : x_offset(NAN), y_offset(NAN), band(NAN),
	    pol_angle(NAN), pol_efficiency(NAN), coupling(Unknown) {}
	std::string Description() const;

	std::string physical_name; /* e.g. D4.A2.3.Y */

	/* Focal plane position in angular units as offsets from boresight */
	double x_offset; /* in az */
	double y_offset; /* and el */

	double band;		/* Standard frequency units */
	double pol_angle;	/* Standard angular units */
	double pol_efficiency;	/* 0-1 */

	enum CouplingType {
		Unknown = 'U',
		Optical = 'O',
		DarkTermination = 'T',
		DarkCrossover = 'X',
		Resistor = 'R',
	};

	CouplingType coupling;  /* Optical coupling type */

	std::string wafer_id;
	std::string pixel_id;
	std::string pixel_type;

	template <class A> void serialize(A &ar, unsigned v);
};

G3_POINTERS(BolometerProperties);
G3_SERIALIZABLE(BolometerProperties, 6);

G3MAP_OF(std::string, BolometerProperties, BolometerPropertiesMap);

#endif

