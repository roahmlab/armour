#ifndef DYNAMICS_CPP
#define DYNAMICS_CPP

#include "Dynamics.h"

KinematicsDynamics::KinematicsDynamics(ConstantAccelerationCurve* traj_input) {
    traj = traj_input;

    // pre-allocate memory
    links = PZsparseArray(NUM_FACTORS * 3, NUM_TIME_STEPS);
    
    // initialize robot properties
    for (int i = 0; i < NUM_JOINTS; i++) {
        trans_matrix(i, 0) = Eigen::MatrixXd::Zero(3, 1);
        trans_matrix(i, 0)(0) = trans[3 * i];
        trans_matrix(i, 0)(1) = trans[3 * i + 1];
        trans_matrix(i, 0)(2) = trans[3 * i + 2];
    }

    trans_matrix(NUM_JOINTS, 0) = Eigen::MatrixXd::Zero(3, 1);
    trans_matrix(NUM_JOINTS, 0)(0) = trans[3 * NUM_JOINTS];
    trans_matrix(NUM_JOINTS, 0)(1) = trans[3 * NUM_JOINTS + 1];
    trans_matrix(NUM_JOINTS, 0)(2) = trans[3 * NUM_JOINTS + 2];

    // define original link PZs
    links = PZsparseArray(NUM_JOINTS, NUM_TIME_STEPS);

    for (int i = 0; i < NUM_JOINTS; i++) {
        PZsparseArray link(3, 1);

        for (int j = 0; j < 3; j++) {
            uint64_t degree[1][NUM_FACTORS * 6] = {0};
            degree[0][NUM_FACTORS * (j + 1)] = 1; // use qde, qdae, qdde for x, y, z generator
            double temp = link_zonotope_generators[i][j];
            link(j, 0) = PZsparse(link_zonotope_center[i][j], &temp, degree, 1);
        }

        links(i, 0) = stack(link);

        for (int j = 1; j < NUM_TIME_STEPS; j++) {
            links(i, j) = links(i, 0);
        }
    }
}

void KinematicsDynamics::fk(uint t_ind) {
    PZsparse FK_R = PZsparse(0, 0, 0); // identity matrix
    PZsparse FK_T(3, 1);
    
    for (int i = 0; i < NUM_JOINTS; i++) {
        PZsparse P(trans_matrix(i, 0));
        
        FK_T = FK_T + FK_R * P;
        FK_R = FK_R * traj->R(i, t_ind);
        
        links(i, t_ind) = FK_R * links(i, t_ind) + FK_T;
    }
}

// void KinematicsDynamics::rnea(uint t_ind,
//                               PZsparseArray& mass_arr,
//                               PZsparseArray& I_arr,
//                               PZsparseArray& u,
//                               bool setGravity) {
//     PZsparse w(3, 1);
//     PZsparse wdot(3, 1);
//     PZsparse w_aux(3, 1);
//     PZsparse linear_acc(3, 1);

//     PZsparseArray F(NUM_JOINTS, 1);
//     PZsparseArray N(NUM_JOINTS, 1);

//     if (setGravity) { // set gravity
//         // directly modify the center of the PZ instance
//         linear_acc.center(2) = gravity;
//     }

//     // RNEA forward recursion
//     for (int i = 0; i < NUM_JOINTS; i++) {
//         // NOTE:
//         // This is just a simplified implementation!!!
//         // We assume all fixed joints are at the end and the revolute joints are consecutive
//         if (axes[i] != 0) { // revolute joints
//             // line 16
//             linear_acc = traj->R_t(i, t_ind) * (linear_acc 
//                                                  + cross(wdot, trans_matrix(i, 0)) 
//                                                  + cross(w, cross(w_aux, trans_matrix(i, 0))));

//             // line 13
//             w = traj->R_t(i, t_ind) * w;
//             w.addOneDimPZ(traj->qd_des(i, t_ind), abs(axes[i]) - 1, 0);

//             // line 14
//             w_aux = traj->R_t(i, t_ind) * w_aux;

//             // line 15
//             wdot = traj->R_t(i, t_ind) * wdot;

//             PZsparse temp(3, 1); // temp = joint_vel(robot_params.q_index(i))*z(:,i)
//             temp.addOneDimPZ(traj->qd_des(i, t_ind), abs(axes[i]) - 1, 0);

//             wdot = wdot + cross(w_aux, temp);

//             wdot.addOneDimPZ(traj->qdda_des(i, t_ind), abs(axes[i]) - 1, 0);

//             // line 14
//             w_aux.addOneDimPZ(traj->qda_des(i, t_ind), abs(axes[i]) - 1, 0);
//         }
//         else { // fixed joints
//             // line 16
//             linear_acc = traj->R_t(i, t_ind) * (linear_acc 
//                                                  + cross(wdot, trans_matrix(i, 0)) 
//                                                  + cross(w, cross(w_aux, trans_matrix(i, 0))));

//             // line 13
//             w = traj->R_t(i, t_ind) * w;

//             // line 14
//             w_aux = traj->R_t(i, t_ind) * w_aux;

//             // line 15
//             wdot = traj->R_t(i, t_ind) * wdot;
//         }

//         // line 23 & 27
//         F(i, 0) = mass_arr(i, 0) * (linear_acc
//                                      + cross(wdot, com_matrix(i, 0))
//                                      + cross(w, cross(w_aux, com_matrix(i, 0))));

//         // line 29
//         N(i, 0) = I_arr(i, 0) * wdot + cross(w_aux, (I_arr(i, 0) * w));
//     }

//     PZsparse f(3, 1);
//     PZsparse n(3, 1);

//     // RNEA reverse recursion
//     for (int i = NUM_JOINTS - 1; i >= 0; i--) {
//         // line 29
//         n = N(i, 0)
//             + traj->R(i + 1, t_ind) * n
//             + cross(com_matrix(i, 0), F(i, 0))
//             + cross(trans_matrix(i + 1, 0), traj->R(i + 1, t_ind) * f);

//         // line 28
//         f = traj->R(i + 1, t_ind) * f + F(i, 0);

//         if (axes[i] != 0) {
//             u(i, t_ind) = n(abs(axes[i]) - 1, 0);

//             u(i, t_ind) = u(i, t_ind) + armature[i] * traj->qdda_des(i, t_ind);

//             u(i, t_ind) = u(i, t_ind) + damping[i] * traj->qd_des(i, t_ind);

//             // friction is directly cut on the torque limits
//         }
//     }
// }

#endif