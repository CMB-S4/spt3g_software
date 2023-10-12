#ifndef _MAPS_POINTING_H
#define _MAPS_POINTING_H

#include <G3Quat.h>

#include <maps/G3SkyMap.h>

// XXX: Define what the following functions do!

G3VectorQuat
get_detector_pointing_quats(double x_offset, double y_offset,
    const G3VectorQuat & trans_quat, MapCoordReference coord_sys);

void get_detector_pointing( 
	double x_offset, double y_offset,
	const G3VectorQuat & trans_quat,
	MapCoordReference coord_sys,
	std::vector<double> & alpha, std::vector<double> & delta);

void get_detector_rotation( 
	double x_offset, double y_offset,
	const G3VectorQuat & trans_quat,
	std::vector<double> & rot);

quat get_transform_quat(double as_0, double ds_0,
			double ae_0, double de_0,
			double as_1, double ds_1,
			double ae_1, double de_1);
quat offsets_to_quat(double x_offset, double y_offset);
quat get_fk5_j2000_to_gal_quat();
void quat_to_ang(quat q, double & alpha, double & delta);
quat ang_to_quat(double alpha, double delta);
quat coord_quat_to_delta_hat(quat q);
quat get_origin_rotator(double alpha, double delta);
double get_rot_ang(quat start_q, quat end_q, quat trans);

double angular_distance(double alpha0, double delta0, double alpha1, double delta1);

#endif

