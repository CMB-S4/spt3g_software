#include <pybindings.h>
#include <G3Constants.h>

PYBINDINGS("core") {
	py::object umod = export_namespace(py::scope(), "G3Units");
#define G3_UNITS_DEF(T) \
	umod.attr(#T) = G3Units::T

	/* Time */
	G3_UNITS_DEF(nanosecond);
	G3_UNITS_DEF(ns);
	G3_UNITS_DEF(nanoseconds);
	G3_UNITS_DEF(microsecond);
	G3_UNITS_DEF(microseconds);
	G3_UNITS_DEF(us);
	G3_UNITS_DEF(millisecond);
	G3_UNITS_DEF(milliseconds);
	G3_UNITS_DEF(ms);
	G3_UNITS_DEF(second);
	G3_UNITS_DEF(seconds);
	G3_UNITS_DEF(sec);
	G3_UNITS_DEF(s);
	G3_UNITS_DEF(minute);
	G3_UNITS_DEF(minutes);
	G3_UNITS_DEF(min);
	G3_UNITS_DEF(hour);
	G3_UNITS_DEF(hours);
	G3_UNITS_DEF(h);
	G3_UNITS_DEF(day);
	G3_UNITS_DEF(days);
	G3_UNITS_DEF(rel);

	/* Frequency */
	G3_UNITS_DEF(Hz);
	G3_UNITS_DEF(hz);
	G3_UNITS_DEF(kHz);
	G3_UNITS_DEF(MHz);
	G3_UNITS_DEF(GHz);

	/* Angle */
	G3_UNITS_DEF(rad);
	G3_UNITS_DEF(radian);
	G3_UNITS_DEF(radians);
	G3_UNITS_DEF(deg);
	G3_UNITS_DEF(degree);
	G3_UNITS_DEF(degrees);
	G3_UNITS_DEF(arcmin);
	G3_UNITS_DEF(arcsec);
	G3_UNITS_DEF(rahour);
	G3_UNITS_DEF(raminute);
	G3_UNITS_DEF(rasecond);
	G3_UNITS_DEF(rahr);

	/* Length */
	G3_UNITS_DEF(meter);
	G3_UNITS_DEF(m);
	G3_UNITS_DEF(meters);
	G3_UNITS_DEF(centimeter);
	G3_UNITS_DEF(cm);
	G3_UNITS_DEF(millimeter);
	G3_UNITS_DEF(mm);
	G3_UNITS_DEF(micron);
	G3_UNITS_DEF(nanometer);
	G3_UNITS_DEF(nm);
	G3_UNITS_DEF(kilometer);
	G3_UNITS_DEF(km);
	G3_UNITS_DEF(au);
	G3_UNITS_DEF(AU);
	G3_UNITS_DEF(parsec);
	G3_UNITS_DEF(pc);
	G3_UNITS_DEF(inch);
	G3_UNITS_DEF(in);
	G3_UNITS_DEF(foot);
	G3_UNITS_DEF(ft);

	/* Power */
	G3_UNITS_DEF(watt);
	G3_UNITS_DEF(W);
	G3_UNITS_DEF(milliwatt);
	G3_UNITS_DEF(mW);
	G3_UNITS_DEF(microwatt);
	G3_UNITS_DEF(uW);
	G3_UNITS_DEF(nanowatt);
	G3_UNITS_DEF(nW);
	G3_UNITS_DEF(picowatt);
	G3_UNITS_DEF(pW);
	G3_UNITS_DEF(attowatt);
	G3_UNITS_DEF(aW);
	G3_UNITS_DEF(horsepower);
	G3_UNITS_DEF(hp);

	/* Flux density */
	G3_UNITS_DEF(jansky);
	G3_UNITS_DEF(Jy);
	G3_UNITS_DEF(millijansky);
	G3_UNITS_DEF(mJy);
	G3_UNITS_DEF(megajansky);
	G3_UNITS_DEF(MJy);

	/* Solid angle */
	G3_UNITS_DEF(sr);
	G3_UNITS_DEF(steradian);
	G3_UNITS_DEF(steradians);
	G3_UNITS_DEF(deg2);
	G3_UNITS_DEF(sqdeg);
	G3_UNITS_DEF(squaredegree);
	G3_UNITS_DEF(squaredegrees);
	G3_UNITS_DEF(arcmin2);
	G3_UNITS_DEF(sqarcmin);
	G3_UNITS_DEF(squarearcmin);

	/* Voltage */
	G3_UNITS_DEF(volt);
	G3_UNITS_DEF(V);
	G3_UNITS_DEF(millivolt);
	G3_UNITS_DEF(mV);
	G3_UNITS_DEF(microvolt);
	G3_UNITS_DEF(uV);

	/* Current */
	G3_UNITS_DEF(ampere);
	G3_UNITS_DEF(amp);
	G3_UNITS_DEF(A);
	G3_UNITS_DEF(milliamp);
	G3_UNITS_DEF(mA);
	G3_UNITS_DEF(microamp);
	G3_UNITS_DEF(uA);
	G3_UNITS_DEF(nanoamp);
	G3_UNITS_DEF(nA);

	/* Resistance */
	G3_UNITS_DEF(ohm);
	G3_UNITS_DEF(mohm);
	G3_UNITS_DEF(milliohm);

	/* Temperature */
	G3_UNITS_DEF(snausage);
	G3_UNITS_DEF(kelvin);
	G3_UNITS_DEF(K);
	G3_UNITS_DEF(millikelvin);
	G3_UNITS_DEF(mK);
	G3_UNITS_DEF(microkelvin);
	G3_UNITS_DEF(uK);
	G3_UNITS_DEF(nanokelvin);
	G3_UNITS_DEF(nK);
	G3_UNITS_DEF(picokelvin);
	G3_UNITS_DEF(pK);
	G3_UNITS_DEF(rankine);
	G3_UNITS_DEF(R);

	/* Pressure */
	G3_UNITS_DEF(bar);
	G3_UNITS_DEF(b);
	G3_UNITS_DEF(millibar);
	G3_UNITS_DEF(mb);
	G3_UNITS_DEF(Pa);
	G3_UNITS_DEF(kPa);

	/* File size */
	G3_UNITS_DEF(byte);
	G3_UNITS_DEF(B);
	G3_UNITS_DEF(bit);
	G3_UNITS_DEF(kilobyte);
	G3_UNITS_DEF(KB);
	G3_UNITS_DEF(megabyte);
	G3_UNITS_DEF(MB);
	G3_UNITS_DEF(gigabyte);
	G3_UNITS_DEF(GB);

	/* Mass */
	G3_UNITS_DEF(gram);
	G3_UNITS_DEF(g);
	G3_UNITS_DEF(kilogram);
	G3_UNITS_DEF(kg);
	G3_UNITS_DEF(milligram);
	G3_UNITS_DEF(mg);

	py::object cmod = export_namespace(py::scope(), "G3Constants");
#define G3_CONSTANTS_DEF(T) \
	cmod.attr(#T) = G3Constants::T

	G3_CONSTANTS_DEF(c);
	G3_CONSTANTS_DEF(h);
	G3_CONSTANTS_DEF(hbar);
	G3_CONSTANTS_DEF(k);
	G3_CONSTANTS_DEF(kb);
	G3_CONSTANTS_DEF(G);
	G3_CONSTANTS_DEF(g0);
	G3_CONSTANTS_DEF(e);
}
