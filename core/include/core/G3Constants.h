#ifndef _G3_CONSTANTS_H
#define _G3_CONSTANTS_H

#include <cmath>
#include <G3Units.h>

namespace G3Constants {
	using namespace G3Units;

	/* Math */
	const double pi = M_PI;
	const double twopi = 2.0 * M_PI;
	const double halfpi = M_PI / 2.0;

	/* Speed of light */
	const double c = 299792458.0 * m / s;

	/* Planck constant */
	const double h = 6.62607015e-34 * W * s * s;
	const double hbar = h / 2.0 / pi;

	/* Boltzmann constant */
	const double k = 1.380649e-23 * W * s / K;
	const double kb = k;

	/* Newtonian gravitational constant */
	const double G = 6.6743e-11 * m * m * m / kg / s / s;

	/* Standard gravitational acceleration */
	const double g0 = 9.80665 * m / s / s;

	/* Elementary charge */
	const double e = 1.602176634e-19 * A * s;
}

#endif
