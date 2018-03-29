#ifndef _G3_UNITS_H
#define _G3_UNITS_H

#include <cmath>

namespace G3Units {
	/* Time: base unit is 10s of nanoseconds */
	const double nanosecond = 0.1;
	const double ns = nanosecond;
	const double nanoseconds = nanosecond;
	const double microsecond = 100;
	const double microseconds = microsecond;
	const double us = microsecond;
	const double millisecond = 1000*microsecond;
	const double milliseconds = millisecond;
	const double ms = millisecond;
	const double second = 1000*millisecond;
	const double seconds = second;
	const double sec = second;
	const double s = second;
	const double minute = 60*second;
	const double minutes = minute;
	const double min = minute;
	const double hour = 60*minute;
	const double hours = hour;
	const double h = hour;
	const double day = 24*hour;
	const double days = day;
	const double rel = 1.2 * second;

	/* Frequency */
	const double Hz = 1/second;
	const double hz = Hz;
	const double kHz = 1/millisecond;
	const double MHz = 1/microsecond;
	const double GHz = 1/nanosecond;

	/* Angles: base unit is radians */
	const double rad = 1;
	const double radian = rad;
	const double radians = rad;
	
	const double deg = M_PI/180.;
	const double degree = deg;
	const double degrees = deg;
	const double arcmin = deg/60.;
	const double arcsec = arcmin/60.;
	const double rahour = 360*deg/24.;
	const double raminute = rahour/60.;
	const double rasecond = raminute/60.;
	const double rahr = rahour;

	/* Length: base unit is millimeter */
	const double meter = 1000;
	const double m = meter;
	const double meters = meter;
	const double centimeter = 1e-2*meter;
	const double cm = centimeter;
	const double millimeter = 1e-3*meter;
	const double mm = millimeter;
	const double micron = 1e-6*meter;
	const double nanometer = 1e-9*meter;
	const double nm = nanometer;
	const double kilometer = 1e3*meter;
	const double km = kilometer;
	const double au = 149597870700*meter;
	const double AU = 149597870700*meter;
	const double parsec = 206264.81*au;
	const double pc = parsec;
	const double inch = 2.54*cm;
	const double in = inch;
	const double foot = 12*inch;
	const double ft = foot;

	/* Power: base unit is watt */
	const double watt = 1;
	const double W = watt;
	const double milliwatt = 1e-3*watt;
	const double mW = milliwatt;
	const double microwatt = 1e-6*watt;
	const double uW = microwatt;
	const double nanowatt = 1e-9*watt;
	const double nW = nanowatt;
	const double picowatt = 1e-12*watt;
	const double pW = picowatt;
	const double attowatt = 1e-18*watt;
	const double aW = attowatt;
	const double horsepower = 745.7*watt;
	const double hp = horsepower;

	/* Voltage: base unit is volts */
	const double volt = 1;
	const double V = volt;
	const double millivolt = 1e-3*volt;
	const double mV = millivolt;
	const double microvolt = 1e-6*volt;
	const double uV = microvolt;

	/* Current: base unit is amps */
	const double ampere = 1;
	const double amp = ampere;
	const double A = amp;
	const double milliamp = 1e-3*amp;
	const double mA = milliamp;
	const double microamp = 1e-6*amp;
	const double uA = microamp;
	const double nanoamp = 1e-9*amp;
	const double nA = nanoamp;

	/* Temperature: base unit is mK */
	const double snausage = 42.0;
	const double kelvin = 1e3;
	const double K = kelvin;
	const double millikelvin = 1e-3*K;
	const double mK = millikelvin;
	const double microkelvin = 1e-6*K;
	const double uK = microkelvin;
	const double nanokelvin = 1e-9*K;
	const double nK = nanokelvin;
	const double picokelvin = 1e-12*K;
	const double pK = picokelvin;
	const double rankine = 9./5. * K;
	const double R = rankine;

	/* Pressure: base unit is bar */
	const double bar = 1;
	const double b = bar;
	const double millibar = b*1e-3;
	const double mb = millibar;

        /* File size: base unit is byte */
        const double byte = 1;
        const double B = byte;
        const double bit = byte / 8.;
        const double kilobyte = 1024*byte;
        const double KB = kilobyte;
        const double megabyte = 1024*kilobyte;
        const double MB = megabyte;
        const double gigabyte = 1024*megabyte;
        const double GB = gigabyte;
}

#endif

