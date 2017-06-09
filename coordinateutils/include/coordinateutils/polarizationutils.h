#ifndef COORDINATEUTILS_POLARIZATION_H
#define COORDINATEUTILS_POLARIZATION_H

/*
 * Polarization rotation utility functions
 */

// Scalar rotations of polarization angles from galactic to equatorial
// coordinates at galactic coordinates (l, b). Polarization angles use
// IAU convention.
double pol_angle_gal_to_eq(double l, double b, double pol_in);

double pol_angle_eq_to_gal(double alpha, double delta, double pol_in);

// As above, except provides transforms of (Q,U) rather polarization angle
// and fraction.
void pol_qu_gal_to_eq(double l, double b, double q_gal, double u_gal,
                      double &q_fk5_out, double &u_fk5_out);

// Same for arrays of polarization angles and positions
void pol_angle_gal_to_eq_arr(double *l, double *b,
    double *pol_in, double *pol_out, size_t size);

void pol_angle_eq_to_gal_arr(double *alpha, double *delta,
    double *pol_in, double *pol_out, size_t size);

// Scalar rotations in polar coordinates to a coordinate system defined in
// equatorial coordinates by a pole at (pole_theta, pole_phi) from standard
// equatorial coordinates relative to Earth. 

double pol_angle_coord_transform(double pole_theta, double pole_phi,
    double pix_theta, double pix_phi, double pol_in);

// The same except that it returns the sin and cos of the angle above
void pol_angle_coord_sin_cos(double pole_theta, double pole_phi,
    double pix_theta, double pix_phi, double *pol_sin, double *pol_cos);

#endif // COORDINATEUTILS_POLARIZATION_H

