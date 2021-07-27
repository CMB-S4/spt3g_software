#ifndef _G3_CONSTANTS_H
#define _G3_CONSTANTS_H

#include <cmath>
#include <G3Units.h>

namespace G3Constants {
	/* Speed of light */
	const double c = 299792458.0 * G3Units::m / G3Units::s;

	/* Planck constant */
	const double h = 6.62607015e-34 * G3Units::W * G3Units::s * G3Units::s;
	const double hbar = h / 2.0 / M_PI;

	/* Boltzmann constant */
	const double k = 1.380649e-23 * G3Units::W * G3Units::s / G3Units::K;
	const double kb = k;

	/* Newtonian gravitational constant */
	const double G = 6.6743e-11 * G3Units::m * G3Units::m * G3Units::m / G3Units::kg / G3Units::s / G3Units::s;

	/* Standard gravitational acceleration */
	const double g0 = 9.80665 * G3Units::m / G3Units::s / G3Units::s;

	/* Elementary charge */
	const double e = 1.602176634e-19 * G3Units::A * G3Units::s;
}

#endif
