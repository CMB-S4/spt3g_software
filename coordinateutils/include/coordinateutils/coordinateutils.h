#ifndef COORDINATEUTILS_H
#define COORDINATEUTILS_H

// General utility functions provided by the coordinate utils library

#include <stdint.h>

// Convert arrays of coordinates between galactic and equatorial coordinates
void coord_gal_to_eq(double *l_in, double *b_in, double *ra_out,
    double *dec_out, size_t size);
void coord_eq_to_gal(double *ra_in, double *dec_in,
    double *l_out, double *b_out, size_t size);


#endif //COORDINATEUTILS_H

