#ifndef DYNAMICS_CPP
#define DYNAMICS_CPP

#include "Dynamics.h"

KinematicsDynamics::KinematicsDynamics(BezierCurve* traj_input) {
    traj = traj_input;

    // pre-allocate memory
    links = PZsparseArray(NUM_FACTORS * 3, NUM_TIME_STEPS);
    mass_nominal_arr = PZsparseArray(NUM_JOINTS, 1);
    mass_uncertain_arr = PZsparseArray(NUM_JOINTS, 1);
    I_nominal_arr = PZsparseArray(NUM_JOINTS, 1);
    I_uncertain_arr = PZsparseArray(NUM_JOINTS, 1);
    u_nom = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    u_nom_int = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    f_c_int = PZsparseArray(1,NUM_TIME_STEPS);
    n_c_int = PZsparseArray(1,NUM_TIME_STEPS);
    f_c_nom = PZsparseArray(1,NUM_TIME_STEPS);
    n_c_nom = PZsparseArray(1,NUM_TIME_STEPS);
    r = PZsparseArray(NUM_FACTORS, 1);
    Mr = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);

    // initialize robot properties
    for (int i = 0; i < NUM_JOINTS; i++) {
        trans_matrix(i, 0) = Eigen::MatrixXd::Zero(3, 1);
        trans_matrix(i, 0)(0) = trans[3 * i];
        trans_matrix(i, 0)(1) = trans[3 * i + 1];
        trans_matrix(i, 0)(2) = trans[3 * i + 2];

        // com_matrix(i, 0) = Eigen::MatrixXd::Zero(3, 1);
        // com_matrix(i, 0)(0) = com[i][0];
        // com_matrix(i, 0)(1) = com[i][1];
        // com_matrix(i, 0)(2) = com[i][2];

        Eigen::MatrixXd mass_matrix(1, 1);
        mass_matrix(0) = mass[i];
        mass_nominal_arr(i) = PZsparse(mass_matrix);
        mass_uncertain_arr(i) = PZsparse(mass_matrix, mass_uncertainty[i]);

        Eigen::Matrix3d inertia_matrix;
        for (int j = 0; j < 9; j++) {
            inertia_matrix(j) = inertia[i * 9 + j]; // This may not be right...
        }
        I_nominal_arr(i) = PZsparse(inertia_matrix);
        I_uncertain_arr(i) = PZsparse(inertia_matrix, inertia_uncertainty[i]);

        if (i < NUM_FACTORS) {
            r(i) = PZsparse(0, Interval(-eps, eps));
        }
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

void KinematicsDynamics::fk(uint s_ind) {
    PZsparse FK_R = PZsparse(0, 0, 0); // identity matrix
    PZsparse FK_T(3, 1);
    int j = 0;

    for (int i = 0; i < NUM_JOINTS; i++) {
        PZsparse P(trans_matrix(i, 0));
        
        FK_T = FK_T + FK_R * P;
        FK_R = FK_R * traj->R(i, s_ind);
        
        links(i, s_ind) = FK_R * links(i, s_ind) + FK_T;
    }
}

void KinematicsDynamics::rnea(uint s_ind,
                              PZsparseArray& mass_arr,
                              PZsparseArray& I_arr,
                              PZsparseArray& u,
                              PZsparseArray& f_c,
                              PZsparseArray& n_c,
                              bool setGravity) {
    PZsparse& cq1 = traj->cos_q_des(0, s_ind);
    PZsparse& cq2 = traj->cos_q_des(1, s_ind);
    PZsparse& cq3 = traj->cos_q_des(2, s_ind);
    PZsparse& cq4 = traj->cos_q_des(3, s_ind);
    PZsparse& cq5 = traj->cos_q_des(4, s_ind);
    PZsparse& cq6 = traj->cos_q_des(5, s_ind);
    PZsparse& cq7 = traj->cos_q_des(6, s_ind);

    PZsparse& sq1 = traj->sin_q_des(0, s_ind);
    PZsparse& sq2 = traj->sin_q_des(1, s_ind);
    PZsparse& sq3 = traj->sin_q_des(2, s_ind);
    PZsparse& sq4 = traj->sin_q_des(3, s_ind);
    PZsparse& sq5 = traj->sin_q_des(4, s_ind);
    PZsparse& sq6 = traj->sin_q_des(5, s_ind);
    PZsparse& sq7 = traj->sin_q_des(6, s_ind);

    PZsparse& qd1 = traj->qd_des(0, s_ind);
    PZsparse& qd2 = traj->qd_des(1, s_ind);
    PZsparse& qd3 = traj->qd_des(2, s_ind);
    PZsparse& qd4 = traj->qd_des(3, s_ind);
    PZsparse& qd5 = traj->qd_des(4, s_ind);
    PZsparse& qd6 = traj->qd_des(5, s_ind);
    PZsparse& qd7 = traj->qd_des(6, s_ind);

    PZsparse& qda1 = traj->qda_des(0, s_ind);
    PZsparse& qda2 = traj->qda_des(1, s_ind);
    PZsparse& qda3 = traj->qda_des(2, s_ind);
    PZsparse& qda4 = traj->qda_des(3, s_ind);
    PZsparse& qda5 = traj->qda_des(4, s_ind);
    PZsparse& qda6 = traj->qda_des(5, s_ind);
    PZsparse& qda7 = traj->qda_des(6, s_ind);

    PZsparse& qdd1 = traj->qdda_des(0, s_ind);
    PZsparse& qdd2 = traj->qdda_des(1, s_ind);
    PZsparse& qdd3 = traj->qdda_des(2, s_ind);
    PZsparse& qdd4 = traj->qdda_des(3, s_ind);
    PZsparse& qdd5 = traj->qdda_des(4, s_ind);
    PZsparse& qdd6 = traj->qdda_des(5, s_ind);
    PZsparse& qdd7 = traj->qdda_des(6, s_ind);

    PZsparse w1(0);
    PZsparse w2(0);
    PZsparse w3(0);
    PZsparse w_aux1(0);
    PZsparse w_aux2(0);
    PZsparse w_aux3(0);
    PZsparse wdot1(0);
    PZsparse wdot2(0);
    PZsparse wdot3(0);
    PZsparse linear_acc1(0);
    PZsparse linear_acc2(0);
    PZsparse linear_acc3(0);

    PZsparse w_new1(0);
    PZsparse w_new2(0);
    PZsparse w_new3(0);
    PZsparse w_aux_new1(0);
    PZsparse w_aux_new2(0);
    PZsparse w_aux_new3(0);
    PZsparse wdot_new1(0);
    PZsparse wdot_new2(0);
    PZsparse wdot_new3(0);
    PZsparse linear_acc_new1(0);
    PZsparse linear_acc_new2(0);
    PZsparse linear_acc_new3(0);

    PZsparse t2(0);
    PZsparse t3(0);
    PZsparse t4(0);
    PZsparse t5(0);
    PZsparse t6(0);
    PZsparse t7(0);
    PZsparse t8(0);
    PZsparse t9(0);
    PZsparse t10(0);
    PZsparse t11(0);
    PZsparse t12(0);
    PZsparse t13(0);
    PZsparse t14(0);
    PZsparse t15(0);
    PZsparse t16(0);
    PZsparse t17(0);
    PZsparse t18(0);
    PZsparse t19(0);
    PZsparse t20(0);
    PZsparse t21(0);
    PZsparse t22(0);

    // joint 1
    w_new3 = qd1;

    w_aux_new3 = qda1;

    wdot_new3 = qdd1;

    linear_acc_new3 = -9.81E+2/1.0E+2;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[0][0]*w_aux2;
    t3 = com[0][0]*w_aux3;
    t4 = com[0][1]*w_aux1;
    t5 = com[0][1]*w_aux3;
    t6 = com[0][2]*w_aux1;
    t7 = com[0][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F1_1 = -mass_arr(0,0)*(-linear_acc1+com[0][1]*wdot3-com[0][2]*wdot2+t11*w2+t12*w3);
    PZsparse F1_2 = mass_arr(0,0)*(linear_acc2+com[0][0]*wdot3-com[0][2]*wdot1+t11*w1-t13*w3);
    PZsparse F1_3 = mass_arr(0,0)*(linear_acc3-com[0][0]*wdot2+com[0][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(0,0)(0,0)*w1;
    t3 = I_arr(0,0)(0,1)*w2;
    t4 = I_arr(0,0)(0,2)*w3;
    t5 = I_arr(0,0)(1,0)*w1;
    t6 = I_arr(0,0)(1,1)*w2;
    t7 = I_arr(0,0)(1,2)*w3;
    t8 = I_arr(0,0)(2,0)*w1;
    t9 = I_arr(0,0)(2,1)*w2;
    t10 = I_arr(0,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N1_1 = I_arr(0,0)(0,0)*wdot1+I_arr(0,0)(0,1)*wdot2+I_arr(0,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N1_2 = I_arr(0,0)(1,0)*wdot1+I_arr(0,0)(1,1)*wdot2+I_arr(0,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N1_3 = I_arr(0,0)(2,0)*wdot1+I_arr(0,0)(2,1)*wdot2+I_arr(0,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 2
    w_new1 = cq2*w1-sq2*w2*3.673205103346574E-6+sq2*w3;
    w_new2 = cq2*w2*(-3.673205103346574E-6)+cq2*w3-sq2*w1;
    w_new3 = qd2-w2-w3*3.673205103346574E-6;

    w_aux_new1 = cq2*w_aux1-sq2*w_aux2*3.673205103346574E-6+sq2*w_aux3;
    w_aux_new2 = cq2*w_aux2*(-3.673205103346574E-6)+cq2*w_aux3-sq2*w_aux1;
    w_aux_new3 = qda2-w_aux2-w_aux3*3.673205103346574E-6;

    wdot_new1 = cq2*wdot1-sq2*wdot2*3.673205103346574E-6+sq2*wdot3-qd2*(cq2*w_aux2*3.673205103346574E-6-cq2*w_aux3+sq2*w_aux1);
    wdot_new2 = cq2*wdot2*(-3.673205103346574E-6)+cq2*wdot3-sq2*wdot1-qd2*(cq2*w_aux1-sq2*w_aux2*3.673205103346574E-6+sq2*w_aux3);
    wdot_new3 = qdd2-wdot2-wdot3*3.673205103346574E-6;

    t2 = -linear_acc1;
    t3 = w_aux3*5.375E-3;
    t4 = wdot1*5.375E-3;
    t5 = wdot3*5.375E-3;
    t6 = w_aux1*w1*5.375E-3;
    t7 = w_aux1*w2*5.375E-3;
    t10 = w_aux2*1.2838E-1;
    t11 = wdot1*1.2838E-1;
    t12 = wdot2*1.2838E-1;
    t13 = w_aux1*w1*1.2838E-1;
    t14 = w_aux1*w3*1.2838E-1;
    t8 = -t6;
    t9 = -t7;
    t15 = t3+t10;
    t16 = t15*w2;
    t17 = t15*w3;
    t20 = t2+t5+t9+t12+t14;
    t18 = -t17;
    t19 = linear_acc3+t4+t13+t16;
    t21 = linear_acc2+t8+t11+t18;
    linear_acc_new1 = -cq2*t20+sq2*t19-sq2*t21*3.673205103346574E-6;
    linear_acc_new2 = cq2*t19-cq2*t21*3.673205103346574E-6+sq2*t20;
    linear_acc_new3 = -linear_acc2-linear_acc3*3.673205103346574E-6-t16*3.673205103346574E-6+t17-wdot1*1.283800197434774E-1+w_aux1*w1*5.374528433928832E-3;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[1][0]*w_aux2;
    t3 = com[1][0]*w_aux3;
    t4 = com[1][1]*w_aux1;
    t5 = com[1][1]*w_aux3;
    t6 = com[1][2]*w_aux1;
    t7 = com[1][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F2_1 = -mass_arr(1,0)*(-linear_acc1+com[1][1]*wdot3-com[1][2]*wdot2+t11*w2+t12*w3);
    PZsparse F2_2 = mass_arr(1,0)*(linear_acc2+com[1][0]*wdot3-com[1][2]*wdot1+t11*w1-t13*w3);
    PZsparse F2_3 = mass_arr(1,0)*(linear_acc3-com[1][0]*wdot2+com[1][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(1,0)(0,0)*w1;
    t3 = I_arr(1,0)(0,1)*w2;
    t4 = I_arr(1,0)(0,2)*w3;
    t5 = I_arr(1,0)(1,0)*w1;
    t6 = I_arr(1,0)(1,1)*w2;
    t7 = I_arr(1,0)(1,2)*w3;
    t8 = I_arr(1,0)(2,0)*w1;
    t9 = I_arr(1,0)(2,1)*w2;
    t10 = I_arr(1,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N2_1 = I_arr(1,0)(0,0)*wdot1+I_arr(1,0)(0,1)*wdot2+I_arr(1,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N2_2 = I_arr(1,0)(1,0)*wdot1+I_arr(1,0)(1,1)*wdot2+I_arr(1,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N2_3 = I_arr(1,0)(2,0)*wdot1+I_arr(1,0)(2,1)*wdot2+I_arr(1,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 3
    w_new1 = cq3*w1-sq3*w2*3.673205103346574E-6-sq3*w3;
    w_new2 = cq3*w2*(-3.673205103346574E-6)-cq3*w3-sq3*w1;
    w_new3 = qd3+w2-w3*3.673205103346574E-6;

    w_aux_new1 = cq3*w_aux1-sq3*w_aux2*3.673205103346574E-6-sq3*w_aux3;
    w_aux_new2 = cq3*w_aux2*(-3.673205103346574E-6)-cq3*w_aux3-sq3*w_aux1;
    w_aux_new3 = qda3+w_aux2-w_aux3*3.673205103346574E-6;

    wdot_new1 = cq3*wdot1-sq3*wdot2*3.673205103346574E-6-sq3*wdot3-qd3*(cq3*w_aux2*3.673205103346574E-6+cq3*w_aux3+sq3*w_aux1);
    wdot_new2 = cq3*wdot2*(-3.673205103346574E-6)-cq3*wdot3-sq3*wdot1+qd3*(-cq3*w_aux1+sq3*w_aux2*3.673205103346574E-6+sq3*w_aux3);
    wdot_new3 = qdd3+wdot2-wdot3*3.673205103346574E-6;

    t2 = -linear_acc1;
    t3 = w_aux2*6.375E-3;
    t4 = wdot1*6.375E-3;
    t5 = wdot2*6.375E-3;
    t6 = w_aux1*w1*6.375E-3;
    t7 = w_aux1*w3*6.375E-3;
    t8 = w_aux3*2.1038E-1;
    t9 = wdot1*2.1038E-1;
    t10 = wdot3*2.1038E-1;
    t11 = w_aux1*w1*2.1038E-1;
    t12 = w_aux1*w2*2.1038E-1;
    t13 = -t8;
    t14 = -t9;
    t15 = -t10;
    t16 = t3+t13;
    t20 = t2+t5+t7+t12+t15;
    t17 = t16*w2;
    t18 = t16*w3;
    t19 = -t18;
    t21 = linear_acc3+t6+t14+t17;
    t22 = linear_acc2+t4+t11+t19;
    linear_acc_new1 = -cq3*t20-sq3*t21-sq3*t22*3.673205103346574E-6;
    linear_acc_new2 = -cq3*t21-cq3*t22*3.673205103346574E-6+sq3*t20;
    linear_acc_new3 = linear_acc2-linear_acc3*3.673205103346574E-6-t17*3.673205103346574E-6+t19+wdot1*6.375772768889642E-3+w_aux1*w1*2.103799765833175E-1;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[2][0]*w_aux2;
    t3 = com[2][0]*w_aux3;
    t4 = com[2][1]*w_aux1;
    t5 = com[2][1]*w_aux3;
    t6 = com[2][2]*w_aux1;
    t7 = com[2][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F3_1 = -mass_arr(2,0)*(-linear_acc1+com[2][1]*wdot3-com[2][2]*wdot2+t11*w2+t12*w3);
    PZsparse F3_2 = mass_arr(2,0)*(linear_acc2+com[2][0]*wdot3-com[2][2]*wdot1+t11*w1-t13*w3);
    PZsparse F3_3 = mass_arr(2,0)*(linear_acc3-com[2][0]*wdot2+com[2][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(2,0)(0,0)*w1;
    t3 = I_arr(2,0)(0,1)*w2;
    t4 = I_arr(2,0)(0,2)*w3;
    t5 = I_arr(2,0)(1,0)*w1;
    t6 = I_arr(2,0)(1,1)*w2;
    t7 = I_arr(2,0)(1,2)*w3;
    t8 = I_arr(2,0)(2,0)*w1;
    t9 = I_arr(2,0)(2,1)*w2;
    t10 = I_arr(2,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N3_1 = I_arr(2,0)(0,0)*wdot1+I_arr(2,0)(0,1)*wdot2+I_arr(2,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N3_2 = I_arr(2,0)(1,0)*wdot1+I_arr(2,0)(1,1)*wdot2+I_arr(2,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N3_3 = I_arr(2,0)(2,0)*wdot1+I_arr(2,0)(2,1)*wdot2+I_arr(2,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 4
    w_new1 = cq4*w1-sq4*w2*3.673205103346573E-6+sq4*w3;
    w_new2 = cq4*w2*(-3.673205103346573E-6)+cq4*w3-sq4*w1;
    w_new3 = qd4-w2-w3*3.673205103346573E-6;

    w_aux_new1 = cq4*w_aux1-sq4*w_aux2*3.673205103346573E-6+sq4*w_aux3;
    w_aux_new2 = cq4*w_aux2*(-3.673205103346573E-6)+cq4*w_aux3-sq4*w_aux1;
    w_aux_new3 = qda4-w_aux2-w_aux3*3.673205103346573E-6;

    wdot_new1 = cq4*wdot1-sq4*wdot2*3.673205103346573E-6+sq4*wdot3-qd4*(cq4*w_aux2*3.673205103346573E-6-cq4*w_aux3+sq4*w_aux1);
    wdot_new2 = cq4*wdot2*(-3.673205103346573E-6)+cq4*wdot3-sq4*wdot1-qd4*(cq4*w_aux1-sq4*w_aux2*3.673205103346573E-6+sq4*w_aux3);
    wdot_new3 = qdd4-wdot2-wdot3*3.673205103346573E-6;

    t2 = -linear_acc1;
    t3 = w_aux3*6.375E-3;
    t4 = wdot1*6.375E-3;
    t5 = wdot3*6.375E-3;
    t6 = w_aux1*w1*6.375E-3;
    t7 = w_aux1*w2*6.375E-3;
    t10 = w_aux2*2.1038E-1;
    t11 = wdot1*2.1038E-1;
    t12 = wdot2*2.1038E-1;
    t13 = w_aux1*w1*2.1038E-1;
    t14 = w_aux1*w3*2.1038E-1;
    t8 = -t6;
    t9 = -t7;
    t15 = t3+t10;
    t16 = t15*w2;
    t17 = t15*w3;
    t20 = t2+t5+t9+t12+t14;
    t18 = -t17;
    t19 = linear_acc3+t4+t13+t16;
    t21 = linear_acc2+t8+t11+t18;
    linear_acc_new1 = -cq4*t20+sq4*t19-sq4*t21*3.673205103346573E-6;
    linear_acc_new2 = cq4*t19-cq4*t21*3.673205103346573E-6+sq4*t20;
    linear_acc_new3 = -linear_acc2-linear_acc3*3.673205103346573E-6-t16*3.673205103346573E-6+t17-wdot1*2.103800234166825E-1+w_aux1*w1*6.374227231110358E-3;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[3][0]*w_aux2;
    t3 = com[3][0]*w_aux3;
    t4 = com[3][1]*w_aux1;
    t5 = com[3][1]*w_aux3;
    t6 = com[3][2]*w_aux1;
    t7 = com[3][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F4_1 = -mass_arr(3,0)*(-linear_acc1+com[3][1]*wdot3-com[3][2]*wdot2+t11*w2+t12*w3);
    PZsparse F4_2 = mass_arr(3,0)*(linear_acc2+com[3][0]*wdot3-com[3][2]*wdot1+t11*w1-t13*w3);
    PZsparse F4_3 = mass_arr(3,0)*(linear_acc3-com[3][0]*wdot2+com[3][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(3,0)(0,0)*w1;
    t3 = I_arr(3,0)(0,1)*w2;
    t4 = I_arr(3,0)(0,2)*w3;
    t5 = I_arr(3,0)(1,0)*w1;
    t6 = I_arr(3,0)(1,1)*w2;
    t7 = I_arr(3,0)(1,2)*w3;
    t8 = I_arr(3,0)(2,0)*w1;
    t9 = I_arr(3,0)(2,1)*w2;
    t10 = I_arr(3,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N4_1 = I_arr(3,0)(0,0)*wdot1+I_arr(3,0)(0,1)*wdot2+I_arr(3,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N4_2 = I_arr(3,0)(1,0)*wdot1+I_arr(3,0)(1,1)*wdot2+I_arr(3,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N4_3 = I_arr(3,0)(2,0)*wdot1+I_arr(3,0)(2,1)*wdot2+I_arr(3,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 5
    w_new1 = cq5*w1-sq5*w2*3.673205103346573E-6-sq5*w3;
    w_new2 = cq5*w2*(-3.673205103346573E-6)-cq5*w3-sq5*w1;
    w_new3 = qd5+w2-w3*3.673205103346573E-6;

    w_aux_new1 = cq5*w_aux1-sq5*w_aux2*3.673205103346573E-6-sq5*w_aux3;
    w_aux_new2 = cq5*w_aux2*(-3.673205103346573E-6)-cq5*w_aux3-sq5*w_aux1;
    w_aux_new3 = qda5+w_aux2-w_aux3*3.673205103346573E-6;

    wdot_new1 = cq5*wdot1-sq5*wdot2*3.673205103346573E-6-sq5*wdot3-qd5*(cq5*w_aux2*3.673205103346573E-6+cq5*w_aux3+sq5*w_aux1);
    wdot_new2 = cq5*wdot2*(-3.673205103346573E-6)-cq5*wdot3-sq5*wdot1+qd5*(-cq5*w_aux1+sq5*w_aux2*3.673205103346573E-6+sq5*w_aux3);
    wdot_new3 = qdd5+wdot2-wdot3*3.673205103346573E-6;

    t2 = -linear_acc1;
    t3 = w_aux2*6.375E-3;
    t4 = wdot1*6.375E-3;
    t5 = wdot2*6.375E-3;
    t6 = w_aux1*w1*6.375E-3;
    t7 = w_aux1*w3*6.375E-3;
    t8 = w_aux3*2.0843E-1;
    t9 = wdot1*2.0843E-1;
    t10 = wdot3*2.0843E-1;
    t11 = w_aux1*w1*2.0843E-1;
    t12 = w_aux1*w2*2.0843E-1;
    t13 = -t8;
    t14 = -t9;
    t15 = -t10;
    t16 = t3+t13;
    t20 = t2+t5+t7+t12+t15;
    t17 = t16*w2;
    t18 = t16*w3;
    t19 = -t18;
    t21 = linear_acc3+t6+t14+t17;
    t22 = linear_acc2+t4+t11+t19;
    linear_acc_new1 = -cq5*t20-sq5*t21-sq5*t22*3.673205103346573E-6;
    linear_acc_new2 = -cq5*t21-cq5*t22*3.673205103346573E-6+sq5*t20;
    linear_acc_new3 = linear_acc2-linear_acc3*3.673205103346573E-6-t17*3.673205103346573E-6+t19+wdot1*6.375765606139691E-3+w_aux1*w1*2.084299765833175E-1;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[4][0]*w_aux2;
    t3 = com[4][0]*w_aux3;
    t4 = com[4][1]*w_aux1;
    t5 = com[4][1]*w_aux3;
    t6 = com[4][2]*w_aux1;
    t7 = com[4][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F5_1 = -mass_arr(4,0)*(-linear_acc1+com[4][1]*wdot3-com[4][2]*wdot2+t11*w2+t12*w3);
    PZsparse F5_2 = mass_arr(4,0)*(linear_acc2+com[4][0]*wdot3-com[4][2]*wdot1+t11*w1-t13*w3);
    PZsparse F5_3 = mass_arr(4,0)*(linear_acc3-com[4][0]*wdot2+com[4][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(4,0)(0,0)*w1;
    t3 = I_arr(4,0)(0,1)*w2;
    t4 = I_arr(4,0)(0,2)*w3;
    t5 = I_arr(4,0)(1,0)*w1;
    t6 = I_arr(4,0)(1,1)*w2;
    t7 = I_arr(4,0)(1,2)*w3;
    t8 = I_arr(4,0)(2,0)*w1;
    t9 = I_arr(4,0)(2,1)*w2;
    t10 = I_arr(4,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N5_1 = I_arr(4,0)(0,0)*wdot1+I_arr(4,0)(0,1)*wdot2+I_arr(4,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N5_2 = I_arr(4,0)(1,0)*wdot1+I_arr(4,0)(1,1)*wdot2+I_arr(4,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N5_3 = I_arr(4,0)(2,0)*wdot1+I_arr(4,0)(2,1)*wdot2+I_arr(4,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 6
    w_new1 = cq6*w1-sq6*w2*3.673205103346572E-6+sq6*w3;
    w_new2 = cq6*w2*(-3.673205103346572E-6)+cq6*w3-sq6*w1;
    w_new3 = qd6-w2-w3*3.673205103346572E-6;

    w_aux_new1 = cq6*w_aux1-sq6*w_aux2*3.673205103346572E-6+sq6*w_aux3;
    w_aux_new2 = cq6*w_aux2*(-3.673205103346572E-6)+cq6*w_aux3-sq6*w_aux1;
    w_aux_new3 = qda6-w_aux2-w_aux3*3.673205103346572E-6;

    wdot_new1 = cq6*wdot1-sq6*wdot2*3.673205103346572E-6+sq6*wdot3-qd6*(cq6*w_aux2*3.673205103346572E-6-cq6*w_aux3+sq6*w_aux1);
    wdot_new2 = cq6*wdot2*(-3.673205103346572E-6)+cq6*wdot3-sq6*wdot1-qd6*(cq6*w_aux1-sq6*w_aux2*3.673205103346572E-6+sq6*w_aux3);
    wdot_new3 = qdd6-wdot2-wdot3*3.673205103346572E-6;

    t2 = -linear_acc1;
    t3 = w_aux2*1.0593E-1;
    t4 = wdot1*1.0593E-1;
    t5 = wdot2*1.0593E-1;
    t6 = w_aux1*w1*1.0593E-1;
    t7 = w_aux1*w3*1.0593E-1;
    t8 = w_aux3*1.750499999999995E-4;
    t9 = wdot1*1.750499999999995E-4;
    t10 = wdot3*1.750499999999995E-4;
    t11 = w_aux1*w1*1.750499999999995E-4;
    t12 = w_aux1*w2*1.750499999999995E-4;
    t13 = -t11;
    t14 = -t12;
    t15 = t3+t8;
    t16 = t15*w2;
    t17 = t15*w3;
    t20 = t2+t5+t7+t10+t14;
    t18 = -t17;
    t19 = linear_acc3+t6+t9+t16;
    t21 = linear_acc2+t4+t13+t18;
    linear_acc_new1 = -cq6*t20+sq6*t19-sq6*t21*3.673205103346572E-6;
    linear_acc_new2 = cq6*t19-cq6*t21*3.673205103346572E-6+sq6*t20;
    linear_acc_new3 = -linear_acc2-linear_acc3*3.673205103346572E-6-t16*3.673205103346572E-6+t17-wdot1*1.059300006429946E-1+w_aux1*w1*1.74660897383402E-4;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[5][0]*w_aux2;
    t3 = com[5][0]*w_aux3;
    t4 = com[5][1]*w_aux1;
    t5 = com[5][1]*w_aux3;
    t6 = com[5][2]*w_aux1;
    t7 = com[5][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F6_1 = -mass_arr(5,0)*(-linear_acc1+com[5][1]*wdot3-com[5][2]*wdot2+t11*w2+t12*w3);
    PZsparse F6_2 = mass_arr(5,0)*(linear_acc2+com[5][0]*wdot3-com[5][2]*wdot1+t11*w1-t13*w3);
    PZsparse F6_3 = mass_arr(5,0)*(linear_acc3-com[5][0]*wdot2+com[5][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(5,0)(0,0)*w1;
    t3 = I_arr(5,0)(0,1)*w2;
    t4 = I_arr(5,0)(0,2)*w3;
    t5 = I_arr(5,0)(1,0)*w1;
    t6 = I_arr(5,0)(1,1)*w2;
    t7 = I_arr(5,0)(1,2)*w3;
    t8 = I_arr(5,0)(2,0)*w1;
    t9 = I_arr(5,0)(2,1)*w2;
    t10 = I_arr(5,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N6_1 = I_arr(5,0)(0,0)*wdot1+I_arr(5,0)(0,1)*wdot2+I_arr(5,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N6_2 = I_arr(5,0)(1,0)*wdot1+I_arr(5,0)(1,1)*wdot2+I_arr(5,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N6_3 = I_arr(5,0)(2,0)*wdot1+I_arr(5,0)(2,1)*wdot2+I_arr(5,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 7
    w_new1 = cq7*w1-sq7*w2*3.673205103346572E-6-sq7*w3;
    w_new2 = cq7*w2*(-3.673205103346572E-6)-cq7*w3-sq7*w1;
    w_new3 = qd7+w2-w3*3.673205103346572E-6;

    w_aux_new1 = cq7*w_aux1-sq7*w_aux2*3.673205103346572E-6-sq7*w_aux3;
    w_aux_new2 = cq7*w_aux2*(-3.673205103346572E-6)-cq7*w_aux3-sq7*w_aux1;
    w_aux_new3 = qda7+w_aux2-w_aux3*3.673205103346572E-6;

    wdot_new1 = cq7*wdot1-sq7*wdot2*3.673205103346572E-6-sq7*wdot3-qd7*(cq7*w_aux2*3.673205103346572E-6+cq7*w_aux3+sq7*w_aux1);
    wdot_new2 = cq7*wdot2*(-3.673205103346572E-6)-cq7*wdot3-sq7*wdot1+qd7*(-cq7*w_aux1+sq7*w_aux2*3.673205103346572E-6+sq7*w_aux3);
    wdot_new3 = qdd7+wdot2-wdot3*3.673205103346572E-6;

    t2 = -linear_acc1;
    t3 = w_aux3*1.0593E-1;
    t4 = wdot1*1.0593E-1;
    t5 = wdot3*1.0593E-1;
    t6 = w_aux1*w1*1.0593E-1;
    t7 = w_aux1*w2*1.0593E-1;
    t11 = w_aux2*1.750499999999995E-4;
    t12 = wdot1*1.750499999999995E-4;
    t13 = wdot2*1.750499999999995E-4;
    t14 = w_aux1*w1*1.750499999999995E-4;
    t15 = w_aux1*w3*1.750499999999995E-4;
    t8 = -t3;
    t9 = -t4;
    t10 = -t5;
    t17 = -w2*(t3-t11);
    t18 = -w3*(t3-t11);
    t19 = w3*(t3-t11);
    t16 = t8+t11;
    t20 = t2+t7+t10+t13+t15;
    t21 = linear_acc3+t9+t14+t17;
    t22 = linear_acc2+t6+t12+t19;
    linear_acc_new1 = -cq7*t20-sq7*t21-sq7*t22*3.673205103346572E-6;
    linear_acc_new2 = -cq7*t21-cq7*t22*3.673205103346572E-6+sq7*t20;
    linear_acc_new3 = linear_acc2-linear_acc3*3.673205103346572E-6+t19+wdot1*1.75439102616597E-4+w_aux1*w1*1.059299993570054E-1+w2*(t3-t11)*3.673205103346572E-6;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[6][0]*w_aux2;
    t3 = com[6][0]*w_aux3;
    t4 = com[6][1]*w_aux1;
    t5 = com[6][1]*w_aux3;
    t6 = com[6][2]*w_aux1;
    t7 = com[6][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F7_1 = -mass_arr(6,0)*(-linear_acc1+com[6][1]*wdot3-com[6][2]*wdot2+t11*w2+t12*w3);
    PZsparse F7_2 = mass_arr(6,0)*(linear_acc2+com[6][0]*wdot3-com[6][2]*wdot1+t11*w1-t13*w3);
    PZsparse F7_3 = mass_arr(6,0)*(linear_acc3-com[6][0]*wdot2+com[6][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(6,0)(0,0)*w1;
    t3 = I_arr(6,0)(0,1)*w2;
    t4 = I_arr(6,0)(0,2)*w3;
    t5 = I_arr(6,0)(1,0)*w1;
    t6 = I_arr(6,0)(1,1)*w2;
    t7 = I_arr(6,0)(1,2)*w3;
    t8 = I_arr(6,0)(2,0)*w1;
    t9 = I_arr(6,0)(2,1)*w2;
    t10 = I_arr(6,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N7_1 = I_arr(6,0)(0,0)*wdot1+I_arr(6,0)(0,1)*wdot2+I_arr(6,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N7_2 = I_arr(6,0)(1,0)*wdot1+I_arr(6,0)(1,1)*wdot2+I_arr(6,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N7_3 = I_arr(6,0)(2,0)*wdot1+I_arr(6,0)(2,1)*wdot2+I_arr(6,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    // joint 8
    w_new1 = w1;
    w_new2 = w2;
    w_new3 = w3;

    w_aux_new1 = w_aux1;
    w_aux_new2 = w_aux2;
    w_aux_new3 = w_aux3;

    wdot_new1 = wdot1;
    wdot_new2 = wdot2;
    wdot_new3 = wdot3;

    t2 = w_aux2/5.0E+2;
    t3 = w_aux3*(5.3E+1/5.0E+2);
    t4 = t2+t3;
    linear_acc_new1 = linear_acc1+wdot2/5.0E+2+wdot3*(5.3E+1/5.0E+2)-w_aux1*w2*(5.3E+1/5.0E+2)+(w_aux1*w3)/5.0E+2;
    linear_acc_new2 = linear_acc2-wdot1/5.0E+2+t4*w3+w_aux1*w1*(5.3E+1/5.0E+2);
    linear_acc_new3 = linear_acc3-wdot1*(5.3E+1/5.0E+2)-t4*w2-(w_aux1*w1)/5.0E+2;

    w1 = w_new1;
    w2 = w_new2;
    w3 = w_new3;
    w_aux1 = w_aux_new1;
    w_aux2 = w_aux_new2;
    w_aux3 = w_aux_new3;
    wdot1 = wdot_new1;
    wdot2 = wdot_new2;
    wdot3 = wdot_new3;
    linear_acc1 = linear_acc_new1;
    linear_acc2 = linear_acc_new2;
    linear_acc3 = linear_acc_new3;

    t2 = com[7][0]*w_aux2;
    t3 = com[7][0]*w_aux3;
    t4 = com[7][1]*w_aux1;
    t5 = com[7][1]*w_aux3;
    t6 = com[7][2]*w_aux1;
    t7 = com[7][2]*w_aux2;
    t8 = -t4;
    t9 = -t6;
    t10 = -t7;
    t11 = t2+t8;
    t12 = t3+t9;
    t13 = t5+t10;
    PZsparse F8_1 = -mass_arr(7,0)*(-linear_acc1+com[7][1]*wdot3-com[7][2]*wdot2+t11*w2+t12*w3);
    PZsparse F8_2 = mass_arr(7,0)*(linear_acc2+com[7][0]*wdot3-com[7][2]*wdot1+t11*w1-t13*w3);
    PZsparse F8_3 = mass_arr(7,0)*(linear_acc3-com[7][0]*wdot2+com[7][1]*wdot1+t12*w1+t13*w2);

    t2 = I_arr(7,0)(0,0)*w1;
    t3 = I_arr(7,0)(0,1)*w2;
    t4 = I_arr(7,0)(0,2)*w3;
    t5 = I_arr(7,0)(1,0)*w1;
    t6 = I_arr(7,0)(1,1)*w2;
    t7 = I_arr(7,0)(1,2)*w3;
    t8 = I_arr(7,0)(2,0)*w1;
    t9 = I_arr(7,0)(2,1)*w2;
    t10 = I_arr(7,0)(2,2)*w3;
    t11 = t2+t3+t4;
    t12 = t5+t6+t7;
    t13 = t8+t9+t10;
    PZsparse N8_1 = I_arr(7,0)(0,0)*wdot1+I_arr(7,0)(0,1)*wdot2+I_arr(7,0)(0,2)*wdot3-t12*w_aux3+t13*w_aux2;
    PZsparse N8_2 = I_arr(7,0)(1,0)*wdot1+I_arr(7,0)(1,1)*wdot2+I_arr(7,0)(1,2)*wdot3+t11*w_aux3-t13*w_aux1;
    PZsparse N8_3 = I_arr(7,0)(2,0)*wdot1+I_arr(7,0)(2,1)*wdot2+I_arr(7,0)(2,2)*wdot3-t11*w_aux2+t12*w_aux1;

    PZsparse f8_1 = F8_1;
    PZsparse f8_2 = F8_2;
    PZsparse f8_3 = F8_3;

    PZsparse n8_1 = N8_1+F8_3*com[7][1]-F8_2*com[7][2];
    PZsparse n8_2 = N8_2-F8_3*com[7][0]+F8_1*com[7][2];
    PZsparse n8_3 = N8_3+F8_2*com[7][0]-F8_1*com[7][1];

    PZsparse f7_1 = F7_1+f8_1;
    PZsparse f7_2 = F7_2+f8_2;
    PZsparse f7_3 = F7_3+f8_3;

    PZsparse n7_1 = N7_1-f8_2/5.0E+2-f8_3*(5.3E+1/5.0E+2)+n8_1+F7_3*com[6][1]-F7_2*com[6][2];
    PZsparse n7_2 = N7_2+f8_1/5.0E+2+n8_2-F7_3*com[6][0]+F7_1*com[6][2];
    PZsparse n7_3 = N7_3+f8_1*(5.3E+1/5.0E+2)+n8_3+F7_2*com[6][0]-F7_1*com[6][1];

    PZsparse f6_1 = F6_1+cq7*f7_1-f7_2*sq7;
    PZsparse f6_2 = F6_2+f7_3-cq7*f7_2*3.673205103346572E-6-f7_1*sq7*3.673205103346572E-6;
    PZsparse f6_3 = F6_3-f7_3*3.673205103346572E-6-cq7*f7_2-f7_1*sq7;

    PZsparse n6_1 = N6_1+f7_3*1.75439102616597E-4+F6_3*com[5][1]-F6_2*com[5][2]+cq7*f7_2*1.059299993570054E-1+cq7*n7_1+f7_1*sq7*1.059299993570054E-1-n7_2*sq7;
    PZsparse n6_2 = N6_2+n7_3-F6_3*com[5][0]+F6_1*com[5][2]-cq7*f7_1*1.750499999999995E-4-cq7*n7_2*3.673205103346572E-6+f7_2*sq7*1.750499999999995E-4-n7_1*sq7*3.673205103346572E-6;
    PZsparse n6_3 = N6_3-n7_3*3.673205103346572E-6+F6_2*com[5][0]-F6_1*com[5][1]+cq7*f7_1*1.0593E-1-cq7*n7_2-f7_2*sq7*1.0593E-1-n7_1*sq7;

    t2 = cq6*f6_2;
    t3 = f6_1*sq6;
    PZsparse f5_1 = F5_1+cq6*f6_1-f6_2*sq6;
    PZsparse f5_2 = F5_2-f6_3-t2*3.673205103346572E-6-t3*3.673205103346572E-6;
    PZsparse f5_3 = F5_3-f6_3*3.673205103346572E-6+t2+t3;

    t2 = cq6*n6_2;
    t3 = n6_1*sq6;
    PZsparse n5_1 = N5_1-f6_3*1.059300006429946E-1+F5_3*com[4][1]-F5_2*com[4][2]+cq6*f6_2*1.74660897383402E-4+cq6*n6_1+f6_1*sq6*1.74660897383402E-4-n6_2*sq6;
    PZsparse n5_2 = N5_2-n6_3-t2*3.673205103346572E-6-t3*3.673205103346572E-6-F5_3*com[4][0]+F5_1*com[4][2]-cq6*f6_1*1.0593E-1+f6_2*sq6*1.0593E-1;
    PZsparse n5_3 = N5_3-n6_3*3.673205103346572E-6+t2+t3+F5_2*com[4][0]-F5_1*com[4][1]-cq6*f6_1*1.750499999999995E-4+f6_2*sq6*1.750499999999995E-4;

    PZsparse f4_1 = F4_1+cq5*f5_1-f5_2*sq5;
    PZsparse f4_2 = F4_2+f5_3-cq5*f5_2*3.673205103346573E-6-f5_1*sq5*3.673205103346573E-6;
    PZsparse f4_3 = F4_3-f5_3*3.673205103346573E-6-cq5*f5_2-f5_1*sq5;

    PZsparse n4_1 = N4_1+f5_3*6.375765606139691E-3+F4_3*com[3][1]-F4_2*com[3][2]+cq5*f5_2*2.084299765833175E-1+cq5*n5_1+f5_1*sq5*2.084299765833175E-1-n5_2*sq5;
    PZsparse n4_2 = N4_2+n5_3-F4_3*com[3][0]+F4_1*com[3][2]-cq5*f5_1*6.375E-3-cq5*n5_2*3.673205103346573E-6+f5_2*sq5*6.375E-3-n5_1*sq5*3.673205103346573E-6;
    PZsparse n4_3 = N4_3-n5_3*3.673205103346573E-6+F4_2*com[3][0]-F4_1*com[3][1]+cq5*f5_1*2.0843E-1-cq5*n5_2-f5_2*sq5*2.0843E-1-n5_1*sq5;

    t2 = cq4*f4_2;
    t3 = f4_1*sq4;
    PZsparse f3_1 = F3_1+cq4*f4_1-f4_2*sq4;
    PZsparse f3_2 = F3_2-f4_3-t2*3.673205103346573E-6-t3*3.673205103346573E-6;
    PZsparse f3_3 = F3_3-f4_3*3.673205103346573E-6+t2+t3;

    t2 = cq4*n4_2;
    t3 = n4_1*sq4;
    PZsparse n3_1 = N3_1-f4_3*2.103800234166825E-1+F3_3*com[2][1]-F3_2*com[2][2]+cq4*f4_2*6.374227231110358E-3+cq4*n4_1+f4_1*sq4*6.374227231110358E-3-n4_2*sq4;
    PZsparse n3_2 = N3_2-n4_3-t2*3.673205103346573E-6-t3*3.673205103346573E-6-F3_3*com[2][0]+F3_1*com[2][2]-cq4*f4_1*2.1038E-1+f4_2*sq4*2.1038E-1;
    PZsparse n3_3 = N3_3-n4_3*3.673205103346573E-6+t2+t3+F3_2*com[2][0]-F3_1*com[2][1]-cq4*f4_1*6.375E-3+f4_2*sq4*6.375E-3;

    PZsparse f2_1 = F2_1+cq3*f3_1-f3_2*sq3;
    PZsparse f2_2 = F2_2+f3_3-cq3*f3_2*3.673205103346574E-6-f3_1*sq3*3.673205103346574E-6;
    PZsparse f2_3 = F2_3-f3_3*3.673205103346574E-6-cq3*f3_2-f3_1*sq3;

    PZsparse n2_1 = N2_1+f3_3*6.375772768889642E-3+F2_3*com[1][1]-F2_2*com[1][2]+cq3*f3_2*2.103799765833175E-1+cq3*n3_1+f3_1*sq3*2.103799765833175E-1-n3_2*sq3;
    PZsparse n2_2 = N2_2+n3_3-F2_3*com[1][0]+F2_1*com[1][2]-cq3*f3_1*6.375E-3-cq3*n3_2*3.673205103346574E-6+f3_2*sq3*6.375E-3-n3_1*sq3*3.673205103346574E-6;
    PZsparse n2_3 = N2_3-n3_3*3.673205103346574E-6+F2_2*com[1][0]-F2_1*com[1][1]+cq3*f3_1*2.1038E-1-cq3*n3_2-f3_2*sq3*2.1038E-1-n3_1*sq3;

    // t2 = cq2*f2_2;
    // t3 = f2_1*sq2;
    // PZsparse f1_1 = F1_1+cq2*f2_1-f2_2*sq2;
    // PZsparse f1_2 = F1_2-f2_3-t2*3.673205103346574E-6-t3*3.673205103346574E-6;
    // PZsparse f1_3 = F1_3-f2_3*3.673205103346574E-6+t2+t3;

    t2 = cq2*n2_2;
    t3 = n2_1*sq2;
    // PZsparse n1_1 = N1_1-f2_3*1.283800197434774E-1+F1_3*com[0][1]-F1_2*com[0][2]+cq2*f2_2*5.374528433928832E-3+cq2*n2_1+f2_1*sq2*5.374528433928832E-3-n2_2*sq2;
    // PZsparse n1_2 = N1_2-n2_3-t2*3.673205103346574E-6-t3*3.673205103346574E-6-F1_3*com[0][0]+F1_1*com[0][2]-cq2*f2_1*1.2838E-1+f2_2*sq2*1.2838E-1;
    PZsparse n1_3 = N1_3-n2_3*3.673205103346574E-6+t2+t3+F1_2*com[0][0]-F1_1*com[0][1]-cq2*f2_1*5.375E-3+f2_2*sq2*5.375E-3;

    u(0,s_ind) = n1_3 + damping[0] * traj->qd_des(0, s_ind) + armature[0] * traj->qdda_des(0, s_ind);
    u(1,s_ind) = n2_3 + damping[1] * traj->qd_des(1, s_ind) + armature[1] * traj->qdda_des(1, s_ind);
    u(2,s_ind) = n3_3 + damping[2] * traj->qd_des(2, s_ind) + armature[2] * traj->qdda_des(2, s_ind);
    u(3,s_ind) = n4_3 + damping[3] * traj->qd_des(3, s_ind) + armature[3] * traj->qdda_des(3, s_ind);
    u(4,s_ind) = n5_3 + damping[4] * traj->qd_des(4, s_ind) + armature[4] * traj->qdda_des(4, s_ind);
    u(5,s_ind) = n6_3 + damping[5] * traj->qd_des(5, s_ind) + armature[5] * traj->qdda_des(5, s_ind);
    u(6,s_ind) = n7_3 + damping[6] * traj->qd_des(6, s_ind) + armature[6] * traj->qdda_des(6, s_ind);

    PZsparseArray f8(3, 1);
    f8(0,0) = f8_1;
    f8(1,0) = f8_2;
    f8(2,0) = f8_3;
    f_c(0,s_ind) = stack(f8); // not sure how to assign these
    PZsparseArray n8(3, 1);
    n8(0,0) = n8_1;
    n8(1,0) = n8_2;
    n8(2,0) = n8_3;
    n_c(0,s_ind) = stack(n8); // not sure how to assign these
}

#endif