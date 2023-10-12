#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "Trajectory.h"

class KinematicsDynamics {
public:
	ConstantAccelerationCurve* traj = nullptr;

	Eigen::Array<Eigen::MatrixXd, NUM_JOINTS + 1, 1> trans_matrix;
	// Eigen::Array<Eigen::MatrixXd, NUM_JOINTS, 1> com_matrix;

	// // inertial paramete PZs
	// // mass_nominal already stored as double array mass in RobotInfo.h
    // PZsparseArray mass_nominal_arr;
	// PZsparseArray mass_uncertain_arr;
    // PZsparseArray I_nominal_arr;
    // PZsparseArray I_uncertain_arr;

	// link PZs
	PZsparseArray links;

	// // nominal torque PZs
    // PZsparseArray u_nom;
    // PZsparseArray u_nom_int;

	KinematicsDynamics() {}

	KinematicsDynamics(ConstantAccelerationCurve* traj_input);

	// generate link PZs through forward kinematics
	void fk(uint t_ind);

	// // generate nominal torque PZs through rnea
	// void rnea(uint t_ind,
	//           PZsparseArray& mass_arr,
	// 		  PZsparseArray& I_arr,
	// 		  PZsparseArray& u,
	// 		  bool setGravity = true);

	// void rnea_nominal(uint t_ind) {
	// 	rnea(t_ind, mass_nominal_arr, I_nominal_arr, u_nom);
	// }

	// void rnea_interval(uint t_ind) {
	// 	rnea(t_ind, mass_uncertain_arr, I_uncertain_arr, u_nom_int);
	// }
};

#endif
