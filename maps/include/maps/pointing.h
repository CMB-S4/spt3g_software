#ifndef _MAPS_POINTING_H
#define _MAPS_POINTING_H

#include <G3Quat.h>

#include <maps/G3SkyMap.h>

// Compute a vector of detector pointing quaternions from a vector of
// boresight transform quaternions
G3VectorQuat
get_detector_pointing_quats(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, MapCoordReference coord_sys);

// Compute a vector of detector pointing pixels from a vector of
// boresight transform quaternions
std::vector<size_t>
get_detector_pointing_pixels(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, G3SkyMapConstPtr skymap);

// Compute vectors of detector pointing angles from a vector of
// boresight transform quaternions
void
get_detector_pointing(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat, MapCoordReference coord_sys,
    std::vector<double> &alpha, std::vector<double> &delta);

// Compute a vector of detector pointing rotation angles from a vector of
// boresight transform quaternions
std::vector<double>
get_detector_rotation(double x_offset, double y_offset,
    const G3VectorQuat &trans_quat);

// Compute a vector quaternion that is the boresight rotated
// by the given x and y offsets.
quat offsets_to_quat(double x_offset, double y_offset);

// Conversion functions between sky coordinates and vector quaternions
void quat_to_ang(quat q, double &alpha, double &delta);
quat ang_to_quat(double alpha, double delta);

// Compute the angular separation between two vector quaternions
double quat_ang_sep(quat q0, quat q1);

// Compute a rotation quaternion that would rotate the boresight vector
// to point in the given sky direction.
quat get_origin_rotator(double alpha, double delta);

// Compute the quaternion for rotating FK5 J2000 coordinates to
// Galactic J2000 coordinates
quat get_fk5_j2000_to_gal_quat();

#endif

