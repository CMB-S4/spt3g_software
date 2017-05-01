#pragma once

double pol_angle_gal_to_eq(double l, double b, double pol_in);

double pol_angle_eq_to_gal(double alpha, double delta, double pol_in);

void pol_angle_gal_to_eq_arr(double *l, double *b, double *pol_in,
                             double *pol_out, size_t size);

void pol_angle_eq_to_gal_arr(double *alpha, double *delta, double *pol_in,
                             double *pol_out, size_t size);

double pol_angle_coord_transform(double pole_theta, double pole_phi,
                                 double pix_theta, double pix_phi, double pol_in);

void pol_angle_coord_sin_cos(double pole_theta, double pole_phi, double pix_theta,
                             double pix_phi, double *pol_sin, double *pol_cos);

void pol_qu_gal_to_eq(double l, double b, double q_gal, double u_gal,
                      double &q_fk5_out, double &u_fk5_out);
