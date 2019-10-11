#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include <coordinateutils/chealpix.h>

/*
  port of get_interp_val from healpix-cxx, because there is no C equivalent.

  to compile:
  gcc -std=c99 get_interp_val.c -o get_interp_val

  ASR, 2016
*/

#ifndef M_PI
#define M_PI		3.14159265358979323846	// pi
#endif
#ifndef M_TWOPI
#define M_TWOPI         6.28318530717958647692  // 2*pi
#endif
#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923	// pi/2
#endif
/*
  Calculate ring parameters for interpolation.

  minfo : map_info structure for all rings
  iring : ring index
  startpix, ringpix, theta, shifted : ring parameters
*/
void get_ring_info(map_info *minfo, long iring, long *startpix,
                   long *ringpix, double *theta, int *shifted) {

  ring_info *ring = minfo->rings + iring;

  if (!ring->init) {
    ring->idx = iring;

    long northring = \
      (iring > (minfo->nring / 2)) ? (minfo->nring - iring) : iring;

    if (northring < minfo->nside) {
      double tmp = northring * northring * minfo->fact2;
      double costheta = 1 - tmp;
      double sintheta = sqrt(tmp * (2 - tmp));
      ring->theta = atan2(sintheta, costheta);
      ring->ringpix = 4 * northring;
      ring->shifted = 1;
      ring->startpix = 2 * northring * (northring - 1);
    } else {
      ring->theta = acos((2 * minfo->nside - northring) * minfo->fact1);
      ring->ringpix = 4 * minfo->nside;
      ring->shifted = (((northring - minfo->nside) & 1) == 0);
      ring->startpix = minfo->ncap + (northring - minfo->nside) * ring->ringpix;
    }

    if (northring != iring) {
      ring->theta = M_PI - ring->theta;
      ring->startpix = minfo->npix - ring->startpix - ring->ringpix;
    }

    ring->init = 1;
  }

  if (theta != NULL)
    *theta = ring->theta;
  if (ringpix != NULL)
    *ringpix = ring->ringpix;
  if (shifted != NULL)
    *shifted = ring->shifted;
  if (startpix != NULL)
    *startpix = ring->startpix;
}

/*
  Allocate and initialize map_info structure.

  nside : map dimension
  shift_phi : if 1, center rings at phi = 0 instead of phi = pi
  populate : if 1, populate all ring parameters
*/
map_info * init_map_info(size_t nside, int nest, int shift_phi, int populate) {
  map_info *minfo = malloc(sizeof(*minfo));
  int iring;
  minfo->nside = nside;
  minfo->npface = minfo->nside * minfo->nside;
  minfo->npix = 12 * minfo->npface;
  minfo->ncap = (minfo->npface - minfo->nside) << 1;
  minfo->nring = 4 * minfo->nside;
  minfo->fact2 = 4. / minfo->npix;
  minfo->fact1 = (minfo->nside << 1) * minfo->fact2;
  minfo->nest = nest;
  minfo->shift_phi = shift_phi;

  minfo->rings = calloc(minfo->nring, sizeof(ring_info));
  if (populate) {
    for (iring = 0; iring < minfo->nring; iring++)
      get_ring_info(minfo, iring, NULL, NULL, NULL, NULL);
  }

  return minfo;
}

/* free map_info structure */
void free_map_info(map_info *minfo) {
  free(minfo->rings);
  free(minfo);
}

/*
  Return the ring index and pixel offset of the given pixel number.
  Assumes ring pixel ordering.

  minfo : initialized map_info structure
  pix : healpix pixel
  iring, ringpix : corresponding ring index and pixel offset within the ring
*/
int get_ring_index(map_info *minfo, long pix, long *iring, long *ringpix) {
  long ncap = minfo->ncap;
  long npix = minfo->npix;
  long nside = minfo->nside;
  long nring = minfo->nring;

  if (minfo->nest)
    nest2ring(nside, pix, &pix);

  if (pix < 0 || pix >= npix)
    return -1;

  if (pix < ncap) /* North Polar cap */
    *iring = (long)(0.5 * (1 + sqrt(1.5 + 2 * pix)));
  else if (pix < (npix - ncap)) /* Equatorial region */
    *iring = (long)((pix - ncap) / nring + nside);
  else /* South Polar cap */
    *iring = nring - (long)(0.5 * (1 + sqrt(2 * (npix - pix) - 0.5)));

  *ringpix = pix - minfo->rings[*iring].startpix;

  if (*ringpix < 0 || *ringpix >= minfo->rings[*iring].ringpix)
    return -1;

  if (minfo->shift_phi)
    *ringpix = (*ringpix + minfo->rings[*iring].ringpix / 2) % minfo->rings[*iring].ringpix;

  return 0;
}

/*
  Return the pixel index of the given ring index and pixel offset.
  Assumes ring pixel ordering.

  minfo : initialized map_info structure
  iring, ringpix : ring index and pixel offset within the ring
  pix : corresponding healpix pixel
 */
int get_pixel_index(map_info *minfo, long iring, long ringpix, long *pix)
{
  if (minfo->shift_phi)
    ringpix = (ringpix + minfo->rings[iring].ringpix / 2) % minfo->rings[iring].ringpix;

  *pix = minfo->rings[iring].startpix + ringpix;

  if (*pix < 0 || *pix > minfo->npix)
    return -1;

  if (minfo->nest)
    ring2nest(minfo->nside, *pix, pix);

  return 0;
}

/*
  Return pixel numbers and weights for bilinear lat/lon interpolation.
  Assumes ring pixel ordering.

  minfo : initialized map_info structure
  theta, phi : interpolation point in radians
  pix, weight: output arrays
*/
int get_interp_weights(map_info *minfo, double theta, double phi,
                       long pix[4], double weight[4]) {
  if (theta < 0 || theta > M_PI)
    return -1;

  int nside = minfo->nside;
  long npix = minfo->npix;
  long nring = minfo->nring;
  double z = cos(theta);
  double az = fabs(z);

  long ir1;
  if (az < 2. / 3.)
    ir1 = nside * (2 - 1.5*z);
  else {
    ir1 = nside * sqrt(3 * (1 - az));
    if (z <= 0)
      ir1 = nring - ir1 - 1;
  }
  long ir2 = ir1 + 1;

  double theta1, theta2, w1, tmp, dphi;
  long sp, nr;
  int shift;
  long i1, i2;

  if (ir1 > 0) {
    get_ring_info(minfo, ir1, &sp, &nr, &theta1, &shift);
    dphi = M_TWOPI / nr;
    tmp = phi / dphi - 0.5 * shift;
    i1 = (tmp < 0) ? tmp - 1 : tmp;
    w1 = (phi - (i1 + 0.5 * shift) * dphi) / dphi;
    if (i1 < 0) i1 += nr;
    i2 = i1 + 1;
    if (i2 >= nr) i2 -= nr;
    pix[0] = sp + i1;
    pix[1] = sp + i2;
    weight[0] = 1 - w1;
    weight[1] = w1;
  }

  if (ir2 < nring) {
    get_ring_info(minfo, ir2, &sp, &nr, &theta2, &shift);
    dphi = M_TWOPI / nr;
    tmp = phi / dphi - 0.5 * shift;
    i1 = (tmp < 0) ? tmp - 1 : tmp;
    w1 = (phi - (i1 + 0.5 * shift) * dphi) / dphi;
    if (i1 < 0) i1 += nr;
    i2 = i1 + 1;
    if (i2 >= nr) i2 -= nr;
    pix[2] = sp + i1;
    pix[3] = sp + i2;
    weight[2] = 1 - w1;
    weight[3] = w1;
  }

  if (ir1 == 0) {
    double wtheta = theta / theta2;
    weight[2] *= wtheta;
    weight[3] *= wtheta;
    double fac = (1 - wtheta) * 0.25;
    weight[0] = fac;
    weight[1] = fac;
    weight[2] += fac;
    weight[3] += fac;
    pix[0] = (pix[2] + 2) & 3;
    pix[1] = (pix[3] + 2) & 3;
  } else if (ir2 == nring) {
    double wtheta = (theta - theta1)/(M_PI - theta1);
    weight[0] *= (1 - wtheta);
    weight[1] *= (1 - wtheta);
    double fac = wtheta * 0.25;
    weight[0] += fac;
    weight[1] += fac;
    weight[2] = fac;
    weight[3] = fac;
    pix[2] = ((pix[0] + 2) & 3) + npix - 4;
    pix[3] = ((pix[1] + 2) & 3) + npix - 4;
  } else {
    double wtheta = (theta - theta1) / (theta2 - theta1);
    weight[0] *= (1 - wtheta);
    weight[1] *= (1 - wtheta);
    weight[2] *= wtheta;
    weight[3] *= wtheta;
  }

  if (minfo->nest) {
    ring2nest(nside, pix[0], pix + 0);
    ring2nest(nside, pix[1], pix + 1);
    ring2nest(nside, pix[2], pix + 2);
    ring2nest(nside, pix[3], pix + 3);
  }

  return 0;
}

/*
  Return interpolated value at the given location on the map.
  Assumes ring pixel ordering.

  minfo : initialized map_info structure
  map : input map, must have nside == minfo->nside
  theta, phi: interpolation point in radians
*/
double get_interp_val(map_info *minfo, double *map,
                      double theta, double phi) {
  long pix[4];
  double weight[4];
  int ii;
  get_interp_weights(minfo, theta, phi, pix, weight);

  double val = 0;
  for (ii = 0; ii < 4; ii++)
    val += map[pix[ii]] * weight[ii];
  return val;
}

/*
  Calculate many interpolated values at the given locations on the map.
  This function allocates and frees the map_info structure for use with
  all the given interpolation points.
  Assumes ring pixel ordering.

  nside : map dimension
  map : map array
  theta, phi : interpolation points in radians
  val : output array of interpolated values
  n : number of interpolation points
*/
void get_interp_valn(int nside, int nest, double *map, double *theta, double *phi,
                     double *val, int n) {
  map_info *minfo = init_map_info(nside, nest, 0, 0);
  int ii;
  for (ii = 0; ii < n; ii++) {
    val[ii] = get_interp_val(minfo, map, theta[ii], phi[ii]);
  }

  free_map_info(minfo);
}

/**
int main(int argc, char *argv[]) {
  int nside = 256;
  long npix = 12 * nside * nside;


  double *map = malloc(npix * sizeof(double));
  for (int ii = 0; ii < npix; ii++) {
    map[ii] = (float) ii;
  }


  double theta[100];
  double phi[100];
  for (int ii = 0; ii < 100; ii++) {
    theta[ii] = M_PI / 100. * ii;
    phi[ii] = M_TWOPI / 100. * ii;
  }


  double val[100];
  get_interp_valn(nside, 0, map, theta, phi, val, 100);

  for (int ii = 0; ii < 100; ii ++) {
    printf("theta: %6f, phi: %6f, val: %6f\n", theta[ii], phi[ii], val[ii]);
  }
}
**/
