#ifndef TRAJECTORY_CPP
#define TRAJECTORY_CPP

#include "Trajectory.h"

BezierCurve::BezierCurve() {
    q0 = Eigen::VectorXd::Zero(NUM_FACTORS);
    Tqd0 = Eigen::VectorXd::Zero(NUM_FACTORS);
    TTqdd0 = Eigen::VectorXd::Zero(NUM_FACTORS);
    Tqd0 = Eigen::VectorXd::Zero(NUM_FACTORS);
    TTqdd0 = Eigen::VectorXd::Zero(NUM_FACTORS);
    ds = 1.0 / NUM_TIME_STEPS;
}

BezierCurve::BezierCurve(const Eigen::VectorXd& q0_inp, 
                         const Eigen::VectorXd& qd0_inp, 
                         const Eigen::VectorXd& qdd0_inp) {
    q0 = q0_inp;
    qd0 = qd0_inp;
    qdd0 = qdd0_inp;   

    Tqd0 = qd0 * DURATION; 
    TTqdd0 = qdd0 * DURATION * DURATION; 

    // cout << TTqdd0.transpose() << endl;

    // pre-allocate memory
    cos_q_des = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    sin_q_des = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    R = PZsparseArray(NUM_JOINTS + 1, NUM_TIME_STEPS);
    R_t = PZsparseArray(NUM_JOINTS, NUM_TIME_STEPS);
    qd_des = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    qda_des = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);
    qdda_des = PZsparseArray(NUM_FACTORS, NUM_TIME_STEPS);

    // initialize the extrema of the k independent part of q_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        q_des_k_indep_extrema_1[i] = (2*Tqd0[i] + TTqdd0[i] + sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2)))/(5*(6*Tqd0[i] + TTqdd0[i]));
        q_des_k_indep_extrema_2[i] = (2*Tqd0[i] + TTqdd0[i] - sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2)))/(5*(6*Tqd0[i] + TTqdd0[i]));
        q_des_k_indep_extremum_1[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], q_des_k_indep_extrema_1[i]);
        q_des_k_indep_extremum_2[i] = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], q_des_k_indep_extrema_2[i]);
    }

    // initialize the extrema of the k independent part of qd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qd_des_k_indep_extrema_1[i] = (18*Tqd0[i] + 4*TTqdd0[i] + sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i]));
        qd_des_k_indep_extrema_2[i] = (18*Tqd0[i] + 4*TTqdd0[i] - sqrt(6*(54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i]));
        qd_des_k_indep_extremum_1[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], qd_des_k_indep_extrema_1[i]);
        qd_des_k_indep_extremum_2[i] = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], qd_des_k_indep_extrema_2[i]);
    }

    // initialize the extrema of the k independent part of qdd_des
    for (int i = 0; i < NUM_FACTORS; i++) {
        qdd_des_k_indep_extrema_1[i] = (32*Tqd0[i] + 6*TTqdd0[i] + sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i]));
        qdd_des_k_indep_extrema_2[i] = (32*Tqd0[i] + 6*TTqdd0[i] - sqrt(2*(152*pow(Tqd0[i],2) + 42*Tqd0[i]*TTqdd0[i] + 3*pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] + TTqdd0[i]));
        qdd_des_k_indep_extremum_1[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], qdd_des_k_indep_extrema_1[i]);
        qdd_des_k_indep_extremum_2[i] = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], qdd_des_k_indep_extrema_2[i]);
    }

    ds = 1.0 / NUM_TIME_STEPS;
}

void BezierCurve::makePolyZono(int s_ind) {
    assert(s_ind < NUM_TIME_STEPS);

    const double s_lb = s_ind * ds;
    const double s_ub = (s_ind + 1) * ds;

    const Interval t_int(s_lb, s_ub);

    for (int i = 0; i < NUM_FACTORS; i++) {
        const double k_range_elt = k_range[i];

        // Part 1: q_des
        double k_dep_coeff_lb = pow(s_lb,3) * (6 * pow(s_lb,2) - 15 * s_lb + 10);
        double k_dep_coeff_ub = pow(s_ub,3) * (6 * pow(s_ub,2) - 15 * s_ub + 10);
        double k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5;
        double k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;
        
        double k_indep_lb = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], s_lb);
        double k_indep_ub = q_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], s_ub);
        if (k_indep_lb > k_indep_ub) {
            swap(k_indep_lb, k_indep_ub);
        }
        if (s_lb < q_des_k_indep_extrema_1[i] && q_des_k_indep_extrema_1[i] < s_ub) {
            k_indep_lb = min(k_indep_lb, q_des_k_indep_extremum_1[i]);
            k_indep_ub = max(k_indep_ub, q_des_k_indep_extremum_1[i]);
        }
        if (s_lb < q_des_k_indep_extrema_2[i] && q_des_k_indep_extrema_2[i] < s_ub) {
            k_indep_lb = min(k_indep_lb, q_des_k_indep_extremum_2[i]);
            k_indep_ub = max(k_indep_ub, q_des_k_indep_extremum_2[i]);
        }
        double k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        double q_des_center = (k_indep_lb + k_indep_ub) * 0.5;
        
        // q_des_k_dep = k_dep_coeff_center * k;
        Interval q_des_radius_int(-k_dep_coeff_radius - k_indep_radius - qe, k_dep_coeff_radius + k_indep_radius + qe);
        
        // q_des_int = q_des_center + q_des_k_dep + q_des_radius_int;

        // first order Taylor expansion
        // Part 1.a: cos(q_des) 
        double cos_q_des_center = cos(q_des_center);
        Interval cos_q_des_radius_int = - q_des_radius_int * sin(q_des_center) 
                                        - 0.5 * cos(q_des_center + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt) + q_des_radius_int) 
                                            * pow(q_des_radius_int + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt), 2);

        cos_q_des_center += getCenter(cos_q_des_radius_int);
        cos_q_des_radius_int = cos_q_des_radius_int - getCenter(cos_q_des_radius_int);
        double cos_q_des_coeff[] = {-k_dep_coeff_center * k_range_elt * sin(q_des_center), getRadius(cos_q_des_radius_int)}; 

        // cos_q_des_int = cos_q_des_center + cos_q_des_coeff[0] * k + cos_q_des_coeff[1] * cosqe;
        uint64_t cos_q_des_degree[2][NUM_FACTORS * 6] = {0};
        cos_q_des_degree[0][i] = 1; // k
        cos_q_des_degree[1][i + NUM_FACTORS * 4] = 1; // cosqe

        cos_q_des(i, s_ind) = PZsparse(cos_q_des_center, cos_q_des_coeff, cos_q_des_degree, 2);

        // Part 1.b: sin(q_des) 
        double sin_q_des_center = sin(q_des_center);
        Interval sin_q_des_radius_int = q_des_radius_int * cos(q_des_center) 
                                        - 0.5 * sin(q_des_center + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt) + q_des_radius_int) 
                                            * pow(q_des_radius_int + k_dep_coeff_center * Interval(-k_range_elt, k_range_elt), 2);

        sin_q_des_center += getCenter(sin_q_des_radius_int);
        sin_q_des_radius_int = sin_q_des_radius_int - getCenter(sin_q_des_radius_int);
        double sin_q_des_coeff[] = {k_dep_coeff_center * k_range_elt * cos(q_des_center), getRadius(sin_q_des_radius_int)};

        // sin_q_des_int = sin_q_des_center + sin_q_des_coeff[0] * k + sin_q_des_coeff[1] * sinqe;
        uint64_t sin_q_des_degree[2][NUM_FACTORS * 6] = {0};
        sin_q_des_degree[0][i] = 1; // k
        sin_q_des_degree[1][i + NUM_FACTORS * 5] = 1; // sinqe

        sin_q_des(i, s_ind) = PZsparse(sin_q_des_center, sin_q_des_coeff, sin_q_des_degree, 2);

        R(i, s_ind) = PZsparse(rots[i * 3], rots[i * 3 + 1], rots[i * 3 + 2]);

        if (axes[i] != 0) {
            R(i, s_ind) = R(i, s_ind) * PZsparse(cos_q_des_center, cos_q_des_coeff, cos_q_des_degree, 2,
                                                 sin_q_des_center, sin_q_des_coeff, sin_q_des_degree, 2,
                                                 axes[i]);
        }

        R_t(i, s_ind) = R(i, s_ind).transpose();
        
        // Part 2: qd_des
        // NOTE:
        // This is just a simplified implementation!!!
        // 30*t^2*(t - 1)^2 in qd_des is just a function with one maxima at t = 0.5
        // So as long as NUM_TIME_STEPS is even number, the following bounding trick holds!
        k_dep_coeff_lb = (30 * pow(s_lb,2) * pow(s_lb - 1,2)) / DURATION;
        k_dep_coeff_ub = (30 * pow(s_ub,2) * pow(s_ub - 1,2)) / DURATION;
        if (k_dep_coeff_ub < k_dep_coeff_lb) { // we are at t >= 0.5, which is a monotonically decreasing region 
            swap(k_dep_coeff_lb, k_dep_coeff_ub);
        }

        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse
        k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse

        k_indep_lb = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], s_lb);
        k_indep_ub = qd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], s_ub);
        if (k_indep_lb > k_indep_ub) {
            swap(k_indep_lb, k_indep_ub);
        }
        if (s_lb < qd_des_k_indep_extrema_1[i] && qd_des_k_indep_extrema_1[i] < s_ub) {
            k_indep_lb = min(k_indep_lb, qd_des_k_indep_extremum_1[i]);
            k_indep_ub = max(k_indep_ub, qd_des_k_indep_extremum_1[i]);
        }
        if (s_lb < qd_des_k_indep_extrema_2[i] && qd_des_k_indep_extrema_2[i] < s_ub) {
            k_indep_lb = min(k_indep_lb, qd_des_k_indep_extremum_2[i]);
            k_indep_ub = max(k_indep_ub, qd_des_k_indep_extremum_2[i]);
        }
        k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        double qd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        double qd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + qde};

        uint64_t qd_des_degree[2][NUM_FACTORS * 6] = {0};
        qd_des_degree[0][i] = 1; // k
        qd_des_degree[1][i + NUM_FACTORS * 1] = 1; // qde

        // qd_des_int = qd_des_center + qd_des_coeff[0] * k + qd_des_coeff[1] * qde;
        qd_des(i, s_ind) = PZsparse(qd_des_center, qd_des_coeff, qd_des_degree, 2);

        double qda_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + qdae};

        uint64_t qda_des_degree[2][NUM_FACTORS * 6] = {0};
        qda_des_degree[0][i] = 1; // k
        qda_des_degree[1][i + NUM_FACTORS * 2] = 1; // qdae

        // qda_des_int = qd_des_center + qda_des_coeff[0] * k + qda_des_coeff[1] * qdae;
        qda_des(i, s_ind) = PZsparse(qd_des_center, qda_des_coeff, qda_des_degree, 2);

        // Part 3: qdd_des
        double temp_lb = (60 * s_lb * (2 * pow(s_lb,2) - 3 * s_lb + 1)) / DURATION / DURATION;
        double temp_ub = (60 * s_ub * (2 * pow(s_ub,2) - 3 * s_ub + 1)) / DURATION / DURATION;
        if (s_ub <= QDD_DES_K_DEP_MAXIMA) { // monotonically increasing region
            k_dep_coeff_lb = temp_lb;
            k_dep_coeff_ub = temp_ub;
        }
        else if (s_lb <= QDD_DES_K_DEP_MAXIMA) { // maxima lives inside
            k_dep_coeff_lb = min(temp_lb, temp_ub);
            k_dep_coeff_ub = (60 * QDD_DES_K_DEP_MAXIMA * (2 * pow(QDD_DES_K_DEP_MAXIMA,2) - 3 * QDD_DES_K_DEP_MAXIMA + 1)) / DURATION / DURATION;
        }
        else if (s_ub <= QDD_DES_K_DEP_MINIMA) { // monotonically decreasing region
            k_dep_coeff_lb = temp_ub;   
            k_dep_coeff_ub = temp_lb;
        }
        else if (s_lb <= QDD_DES_K_DEP_MINIMA) { // minima lives inside
            k_dep_coeff_lb = (60 * QDD_DES_K_DEP_MINIMA * (2 * pow(QDD_DES_K_DEP_MINIMA,2) - 3 * QDD_DES_K_DEP_MINIMA + 1)) / DURATION / DURATION;
            k_dep_coeff_ub = max(temp_lb, temp_ub);
        }
        else { // monotonically increasing region
            k_dep_coeff_lb = temp_lb;
            k_dep_coeff_ub = temp_ub;
        }
        
        k_dep_coeff_center = (k_dep_coeff_ub + k_dep_coeff_lb) * 0.5 * k_range_elt; // Have to scale to [-1,1] in order to fit in PZsparse
        k_dep_coeff_radius = (k_dep_coeff_ub - k_dep_coeff_lb) * 0.5 * k_range_elt;

        k_indep_lb = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], s_lb);
        k_indep_ub = qdd_des_k_indep(q0[i], Tqd0[i], TTqdd0[i], s_ub);
        if (k_indep_lb > k_indep_ub) {
            swap(k_indep_lb, k_indep_ub);
        }
        if (s_lb < qdd_des_k_indep_extrema_1[i] && qdd_des_k_indep_extrema_1[i] < s_ub) {
            k_indep_lb = min(k_indep_lb, qdd_des_k_indep_extremum_1[i]);
            k_indep_ub = max(k_indep_ub, qdd_des_k_indep_extremum_1[i]);
        }
        if (s_lb < qdd_des_k_indep_extrema_2[i] && qdd_des_k_indep_extrema_2[i] < s_ub) {
            k_indep_lb = min(k_indep_lb, qdd_des_k_indep_extremum_2[i]);
            k_indep_ub = max(k_indep_ub, qdd_des_k_indep_extremum_2[i]);
        }
        k_indep_radius = (k_indep_ub - k_indep_lb) * 0.5;
        double qdd_des_center = (k_indep_lb + k_indep_ub) * 0.5;

        double qdd_des_coeff[] = {k_dep_coeff_center, k_dep_coeff_radius + k_indep_radius + qddae};

        uint64_t qdd_des_degree[2][NUM_FACTORS * 6] = {0};
        qdd_des_degree[0][i] = 1; // k
        qdd_des_degree[1][i + NUM_FACTORS * 3] = 1; // qddae

        // qdd_des_int = qdd_des_center + qdd_des_coeff[0] * k + qdd_des_coeff[1] * qdde;
        qdda_des(i, s_ind) = PZsparse(qdd_des_center, qdd_des_coeff, qdd_des_degree, 2);
    }

    // assume all fixed joints are at the end of the kinematics chain
    for (int i = NUM_FACTORS; i < NUM_JOINTS; i++) {
        R(i, s_ind) = PZsparse(rots[i * 3], rots[i * 3 + 1], rots[i * 3 + 2]);
        R_t(i, s_ind) = R(i, s_ind).transpose();
    }

    R(NUM_JOINTS, s_ind) = PZsparse(0, 0, 0);
}

void BezierCurve::returnJointPositionExtremum(double* extremum, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_range[i] * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (2*Tqd0[i] + TTqdd0[i] + sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (2*Tqd0[i] + TTqdd0[i] - sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        double extremum2 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        double extremum3 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        double extremum4 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        double minPosition = min(extremum1, extremum4);
        double maxPosition = max(extremum1, extremum4);
        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            minPosition = min(minPosition, extremum2);
            maxPosition = max(maxPosition, extremum2);
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            minPosition = min(minPosition, extremum3);
            maxPosition = max(maxPosition, extremum3);
        }

        extremum[i              ] = minPosition;
        extremum[i + NUM_FACTORS] = maxPosition;
    }
}

void BezierCurve::returnJointPositionExtremumGradient(double* extremumGradient, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_range[i] * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (2*Tqd0[i] + TTqdd0[i] + sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (2*Tqd0[i] + TTqdd0[i] - sqrt(64*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] - 120*k_actual*Tqd0[i] + pow(TTqdd0[i],2))) / (5*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        double extremum2 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        double extremum3 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        double extremum4 = q_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        double minPosition;
        int minId;
        double maxPosition;
        int maxId;

        if (extremum1 < extremum4) {
            minPosition = extremum1;
            minId = 1;

            maxPosition = extremum4;
            maxId = 4;
        }
        else {
            minPosition = extremum4;
            minId = 4;

            maxPosition = extremum1;
            maxId = 1;
        }

        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum2 < minPosition) {
                minPosition = extremum2;
                minId = 2;
            }
            if (maxPosition < extremum2) {
                maxPosition = extremum2;
                maxId = 2;
            }
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum3 < minPosition) {
                minPosition = extremum3;
                minId = 3;
            }
            if (maxPosition < extremum3) {
                maxPosition = extremum3;
                maxId = 3;
            }
        }

        double minPositionGradient;
        double maxPositionGradient;

        switch (minId) {
            case 1: // t = 0
                minPositionGradient = 0.0;
                break;
            case 2: // t = extrema2
                minPositionGradient = q_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                minPositionGradient = q_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                minPositionGradient = 1.0;
                break;
            default:
                break;
        }

        switch (maxId) {
            case 1: // t = 0
                maxPositionGradient = 0.0;
                break;
            case 2: // t = extrema2
                maxPositionGradient = q_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                maxPositionGradient = q_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                maxPositionGradient = 1.0;
                break;
            default:
                break;
        }

        for (int j = 0; j < NUM_FACTORS; j++) {
            if (i == j) {
                extremumGradient[(i              ) * NUM_FACTORS + j] = minPositionGradient * k_range[i];
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = maxPositionGradient * k_range[i];
            }
            else {
                extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
            }
        }
    }
}

void BezierCurve::returnJointVelocityExtremum(double* extremum, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_range[i] * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] + sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] - sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        double extremum2 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        double extremum3 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        double extremum4 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        double minVelocity = min(extremum1, extremum4);
        double maxVelocity = max(extremum1, extremum4);
        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            minVelocity = min(minVelocity, extremum2); 
            maxVelocity = max(maxVelocity, extremum2);
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            minVelocity = min(minVelocity, extremum3);
            maxVelocity = max(maxVelocity, extremum3);
        }

        extremum[i              ] = minVelocity / DURATION;
        extremum[i + NUM_FACTORS] = maxVelocity / DURATION;
    }
}

void BezierCurve::returnJointVelocityExtremumGradient(double* extremumGradient, const double* k) const {
    for (int i = 0; i < NUM_FACTORS; i++) {
        // k[i] range is [-1,1] since it is defined for PZ, the following is the actual k
        double k_actual = k_range[i] * k[i];

        // list all possible extremas
        double extrema1 = 0;
        double extrema2 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] + sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema3 = (18*Tqd0[i] - 30*k_actual + 4*TTqdd0[i] - sqrt(6*(150*pow(k_actual,2) - 180*k_actual*Tqd0[i] - 20*k_actual*TTqdd0[i] + 54*pow(Tqd0[i],2) + 14*Tqd0[i]*TTqdd0[i] + pow(TTqdd0[i],2))))/(10*(6*Tqd0[i] - 12*k_actual + TTqdd0[i]));
        double extrema4 = 1;

        // get extremums of all extremas
        double extremum1 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema1);
        double extremum2 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema2);
        double extremum3 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema3);
        double extremum4 = qd_des_func(q0[i], Tqd0[i], TTqdd0[i], k_actual, extrema4);

        // find the min and max values
        double minVelocity;
        int minId;
        double maxVelocity;
        int maxId;

        if (extremum1 < extremum4) {
            minVelocity = extremum1;
            minId = 1;

            maxVelocity = extremum4;
            maxId = 4;
        }
        else {
            minVelocity = extremum4;
            minId = 4;

            maxVelocity = extremum1;
            maxId = 1;
        }

        if (0 <= extrema2 && extrema2 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum2 < minVelocity) {
                minVelocity = extremum2;
                minId = 2;
            }
            if (maxVelocity < extremum2) {
                maxVelocity = extremum2;
                maxId = 2;
            }
        }
        if (0 <= extrema3 && extrema3 <= 1) { // check if this extrema is inside the time range [0,1]
            if (extremum3 < minVelocity) {
                minVelocity = extremum3;
                minId = 3;
            }
            if (maxVelocity < extremum3) {
                maxVelocity = extremum3;
                maxId = 3;
            }
        }

        double minVelocityGradient;
        double maxVelocityGradient;

        switch (minId) {
            case 1: // t = 0
                minVelocityGradient = 0.0;
                break;
            case 2: // t = extrema2
                minVelocityGradient = qd_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                minVelocityGradient = qd_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                minVelocityGradient = 1.0;
                break;
            default:
                break;
        }

        switch (maxId) {
            case 1: // t = 0
                maxVelocityGradient = 0.0;
                break;
            case 2: // t = extrema2
                maxVelocityGradient = qd_des_extrema2_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 3: // t = extrema3
                maxVelocityGradient = qd_des_extrema3_k_derivative(q0[i], Tqd0[i], TTqdd0[i], k_actual);
                break;
            case 4: // t = 1
                maxVelocityGradient = 1.0;
                break;
            default:
                break;
        }

        for (int j = 0; j < NUM_FACTORS; j++) {
            if (i == j) {
                extremumGradient[(i              ) * NUM_FACTORS + j] = minVelocityGradient * k_range[i] / DURATION;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = maxVelocityGradient * k_range[i] / DURATION;
            }
            else {
                extremumGradient[(i              ) * NUM_FACTORS + j] = 0.0;
                extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + j] = 0.0;
            }
        }
    }
}

double q_des_func(double q0, double Tqd0, double TTqdd0, double k, double t) {
    double B0 = -pow(t - 1,5);
    double B1 = 5*t*pow(t - 1,4);
    double B2 = -10*pow(t,2)*pow(t - 1,3);
    double B3 = 10*pow(t,3)*pow(t - 1,2);
    double B4 = -5*pow(t,4)*(t - 1);
    double B5 = pow(t,5);
    double beta0 = q0;
    double beta1 = q0 + Tqd0/5;
    double beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    double beta3 = q0 + k;
    double beta4 = q0 + k;
    double beta5 = q0 + k;
    return B0 * beta0 + B1 * beta1 + B2 * beta2 + B3 * beta3 + B4 * beta4 + B5 * beta5;
}

double qd_des_func(double q0, double Tqd0, double TTqdd0, double k, double t) {
    double dB0 = pow(t-1.0,4.0)*-5.0;
    double dB1 = t*pow(t-1.0,3.0)*2.0E+1+pow(t-1.0,4.0)*5.0;
    double dB2 = t*pow(t-1.0,3.0)*-2.0E+1-(t*t)*pow(t-1.0,2.0)*3.0E+1;
    double dB3 = pow(t,3.0)*(t*2.0-2.0)*1.0E+1+(t*t)*pow(t-1.0,2.0)*3.0E+1;
    double dB4 = pow(t,3.0)*(t-1.0)*-2.0E+1-pow(t,4.0)*5.0;
    double dB5 = pow(t,4.0)*5.0;
    double beta0 = q0;
    double beta1 = q0 + Tqd0/5;
    double beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    double beta3 = q0 + k;
    double beta4 = q0 + k;
    double beta5 = q0 + k;
    return dB0 * beta0 + dB1 * beta1 + dB2 * beta2 + dB3 * beta3 + dB4 * beta4 + dB5 * beta5;
}

double qdd_des_func(double q0, double Tqd0, double TTqdd0, double k, double t) {
    double t2 = t*2.0;
    double t3 = t*t;
    double t4 = t*t*t;
    double t5 = t-1.0;
    double t6 = t2-2.0;
    double t7 = t4*2.0E+1;
    double t8 = t5*t5;
    double t9 = t5*t5*t5;
    double t10 = t9*2.0E+1;
    double t11 = t*t8*6.0E+1;
    double t12 = -t10;
    double ddB0 = t12;
    double ddB1 = t9*4.0E+1+t11;
    double ddB2 = t12-t*t8*1.2E+2-t3*t6*3.0E+1;
    double ddB3 = t7+t11+t3*t6*6.0E+1;
    double ddB4 = t4*-4.0E+1-t3*t5*6.0E+1;
    double ddB5 = t7;
    double beta0 = q0;
    double beta1 = q0 + Tqd0/5;
    double beta2 = q0 + (2*Tqd0)/5 + TTqdd0/20;
    double beta3 = q0 + k;
    double beta4 = q0 + k;
    double beta5 = q0 + k;
    return ddB0 * beta0 + ddB1 * beta1 + ddB2 * beta2 + ddB3 * beta3 + ddB4 * beta4 + ddB5 * beta5;
}

double q_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = Tqd0*2.0;
    double t4 = Tqd0*6.0;
    double t5 = Tqd0*Tqd0;
    double t6 = TTqdd0*TTqdd0;
    double t7 = k*1.2E+1;
    double t8 = Tqd0*TTqdd0*1.4E+1;
    double t10 = Tqd0/5.0;
    double t11 = Tqd0*(2.0/5.0);
    double t12 = k*Tqd0*1.2E+2;
    double t13 = TTqdd0/2.0E+1;
    double t9 = -t7;
    double t14 = t5*6.4E+1;
    double t15 = -t12;
    double t16 = q0+t10;
    double t18 = q0+t11+t13;
    double t17 = TTqdd0+t4+t9;
    double t24 = t6+t8+t14+t15;
    double t19 = 1.0/t17;
    double t25 = sqrt(t24);
    double t20 = t19*t19;
    double t21 = t19*t19*t19;
    double t23 = t19*t19*t19*t19*t19;
    double t26 = 1.0/t25;
    double t27 = TTqdd0+t3+t25;
    double t22 = t20*t20;
    double t28 = t27*t27;
    double t29 = t27*t27*t27;
    double t31 = t27*t27*t27*t27*t27;
    double t32 = Tqd0*t19*t26*1.2E+1;
    double t34 = (t19*t27)/5.0;
    double t35 = t20*t27*(1.2E+1/5.0);
    double t30 = t28*t28;
    double t33 = -t32;
    double t36 = t34-1.0;
    double t37 = t36*t36;
    double t38 = t36*t36*t36;
    double t40 = t33+t35;
    double t39 = t37*t37;
    return (t23*t31)/3.125E+3+t2*(t20*t20*t20)*t31*(1.2E+1/6.25E+2)+t21*t29*t37*(2.0/2.5E+1)-(t22*t30*t36)/1.25E+2+q0*t39*(t32-t35)*5.0+t2*t22*t29*t37*(7.2E+1/2.5E+1)-t2*t23*t30*t36*(4.8E+1/1.25E+2)+t16*t20*t27*t39*1.2E+1-t18*t21*t28*t38*(4.8E+1/5.0)+(t2*t22*t30*(t32-t35))/1.25E+2-Tqd0*t2*t23*t26*t30*(1.2E+1/1.25E+2)-Tqd0*t16*t19*t26*t39*6.0E+1-t2*t21*t29*t36*(t32-t35)*(4.0/2.5E+1)-t16*t19*t27*t38*(t32-t35)*4.0+t18*t20*t28*t37*(t32-t35)*(6.0/5.0)-Tqd0*t2*t21*t26*t28*t37*(7.2E+1/5.0)+Tqd0*t2*t22*t26*t29*t36*(4.8E+1/2.5E+1)+Tqd0*t18*t20*t26*t27*t38*4.8E+1;
}

double q_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = Tqd0*2.0;
    double t4 = Tqd0*6.0;
    double t5 = Tqd0*Tqd0;
    double t6 = TTqdd0*TTqdd0;
    double t7 = k*1.2E+1;
    double t8 = Tqd0*TTqdd0*1.4E+1;
    double t10 = Tqd0/5.0;
    double t11 = Tqd0*(2.0/5.0);
    double t12 = k*Tqd0*1.2E+2;
    double t13 = TTqdd0/2.0E+1;
    double t9 = -t7;
    double t14 = t5*6.4E+1;
    double t15 = -t12;
    double t16 = q0+t10;
    double t18 = q0+t11+t13;
    double t17 = TTqdd0+t4+t9;
    double t24 = t6+t8+t14+t15;
    double t19 = 1.0/t17;
    double t25 = sqrt(t24);
    double t20 = t19*t19;
    double t21 = t19*t19*t19;
    double t23 = t19*t19*t19*t19*t19;
    double t26 = 1.0/t25;
    double t27 = -t25;
    double t22 = t20*t20;
    double t28 = TTqdd0+t3+t27;
    double t33 = Tqd0*t19*t26*1.2E+1;
    double t29 = t28*t28;
    double t30 = t28*t28*t28;
    double t32 = t28*t28*t28*t28*t28;
    double t34 = (t19*t28)/5.0;
    double t35 = t20*t28*(1.2E+1/5.0);
    double t31 = t29*t29;
    double t36 = t34-1.0;
    double t40 = t33+t35;
    double t37 = t36*t36;
    double t38 = t36*t36*t36;
    double t39 = t37*t37;
    return (t23*t32)/3.125E+3+t2*(t20*t20*t20)*t32*(1.2E+1/6.25E+2)-q0*t39*t40*5.0+t21*t30*t37*(2.0/2.5E+1)-(t22*t31*t36)/1.25E+2+t2*t22*t30*t37*(7.2E+1/2.5E+1)-t2*t23*t31*t36*(4.8E+1/1.25E+2)-(t2*t22*t31*t40)/1.25E+2+t16*t20*t28*t39*1.2E+1-t18*t21*t29*t38*(4.8E+1/5.0)+Tqd0*t2*t23*t26*t31*(1.2E+1/1.25E+2)+Tqd0*t16*t19*t26*t39*6.0E+1+t2*t21*t30*t36*t40*(4.0/2.5E+1)+t16*t19*t28*t38*t40*4.0-t18*t20*t29*t37*t40*(6.0/5.0)+Tqd0*t2*t21*t26*t29*t37*(7.2E+1/5.0)-Tqd0*t2*t22*t26*t30*t36*(4.8E+1/2.5E+1)-Tqd0*t18*t20*t26*t28*t38*4.8E+1;
}

double qd_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = k*k;
    double t4 = Tqd0*6.0;
    double t5 = TTqdd0*4.0;
    double t6 = Tqd0*Tqd0;
    double t7 = TTqdd0*TTqdd0;
    double t8 = k*1.2E+1;
    double t9 = k*3.0E+1;
    double t10 = Tqd0*1.8E+1;
    double t11 = TTqdd0*2.0E+1;
    double t13 = Tqd0*TTqdd0*1.4E+1;
    double t14 = sqrt(6.0);
    double t17 = k*3.0E+2;
    double t18 = Tqd0*1.8E+2;
    double t21 = k*TTqdd0*-2.0E+1;
    double t24 = k*Tqd0*-1.8E+2;
    double t12 = k*t11;
    double t15 = -t8;
    double t16 = -t9;
    double t19 = t6*5.4E+1;
    double t20 = k*t18;
    double t22 = -t17;
    double t23 = t3*1.5E+2;
    double t25 = TTqdd0+t4+t15;
    double t31 = t11+t18+t22;
    double t32 = t7+t13+t19+t21+t23+t24;
    double t26 = 1.0/t25;
    double t33 = sqrt(t32);
    double t27 = t26*t26;
    double t28 = t26*t26*t26;
    double t30 = t26*t26*t26*t26*t26;
    double t34 = 1.0/t33;
    double t35 = t14*t33;
    double t29 = t27*t27;
    double t36 = t5+t10+t16+t35;
    double t40 = (t14*t31*t34)/2.0;
    double t37 = t36*t36;
    double t38 = t36*t36*t36;
    double t41 = t40+3.0E+1;
    double t42 = (t26*t36)/5.0;
    double t43 = t27*t36*(6.0/5.0);
    double t44 = (t26*t36)/1.0E+1;
    double t39 = t37*t37;
    double t45 = t42-2.0;
    double t46 = t44-1.0;
    double t49 = (t26*t41)/1.0E+1;
    double t47 = t46*t46;
    double t48 = t46*t46*t46;
    double t50 = -t49;
    double t51 = t27*t36*t48*2.4E+1;
    double t52 = t28*t37*t47*(3.6E+1/5.0);
    double t53 = t43+t50;
    double t54 = t26*t41*t48*2.0;
    double t56 = t27*t36*t41*t47*(3.0/5.0);
    double t55 = -t54;
    double t57 = -t56;
    double t58 = t26*t36*t47*t53*6.0;
    double t59 = t27*t37*t46*t53*(3.0/5.0);
    return (q0+Tqd0/5.0)*(t51+t55+t58+t48*t53*2.0E+1)+t2*(t52+t57+t59+(t28*t38*(t27*t36*(1.2E+1/5.0)-(t26*t41)/5.0))/1.0E+2+t29*t38*t45*(9.0/2.5E+1)-t28*t37*t41*t45*(3.0/1.0E+2))-t2*(t30*t39*(3.0/1.25E+2)-(t29*t38*t41)/5.0E+2+t29*t38*t46*(1.8E+1/2.5E+1)+(t28*t38*t53)/5.0E+1-t28*t37*t41*t46*(3.0/5.0E+1))-(q0+Tqd0*(2.0/5.0)+TTqdd0/2.0E+1)*(t51+t52+t55+t57+t58+t59)-q0*t48*t53*2.0E+1+t2*t30*t39*(3.0/1.25E+2)+t27*t37*t47*(3.0/1.0E+1)+(t28*t38*t45)/1.0E+2-(t28*t38*t46)/5.0E+1-(t2*t29*t38*t41)/5.0E+2;
}

double qd_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k) {
    double t2 = k+q0;
    double t3 = k*k;
    double t4 = Tqd0*6.0;
    double t5 = TTqdd0*4.0;
    double t6 = Tqd0*Tqd0;
    double t7 = TTqdd0*TTqdd0;
    double t8 = k*1.2E+1;
    double t9 = k*3.0E+1;
    double t10 = Tqd0*1.8E+1;
    double t12 = TTqdd0*2.0E+1;
    double t14 = Tqd0*TTqdd0*1.4E+1;
    double t15 = sqrt(6.0);
    double t17 = k*3.0E+2;
    double t19 = Tqd0*1.8E+2;
    double t22 = k*TTqdd0*-2.0E+1;
    double t25 = k*Tqd0*-1.8E+2;
    double t11 = -t5;
    double t13 = k*t12;
    double t16 = -t8;
    double t18 = -t10;
    double t20 = t6*5.4E+1;
    double t21 = k*t19;
    double t23 = -t17;
    double t24 = t3*1.5E+2;
    double t26 = TTqdd0+t4+t16;
    double t32 = t12+t19+t23;
    double t33 = t7+t14+t20+t22+t24+t25;
    double t27 = 1.0/t26;
    double t34 = sqrt(t33);
    double t28 = t27*t27;
    double t29 = t27*t27*t27;
    double t31 = t27*t27*t27*t27*t27;
    double t35 = 1.0/t34;
    double t36 = t15*t34;
    double t30 = t28*t28;
    double t37 = t9+t11+t18+t36;
    double t38 = pow(t5-t9+t10-t36,2.0);
    double t39 = -pow(t5-t9+t10-t36,3.0);
    double t41 = (t15*t32*t35)/2.0;
    double t43 = t27*(t5-t9+t10-t36)*(-1.0/5.0);
    double t44 = t28*(t5-t9+t10-t36)*(-6.0/5.0);
    double t45 = t27*(t5-t9+t10-t36)*(-1.0/1.0E+1);
    double t48 = pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,2.0);
    double t49 = -pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0);
    double t52 = t28*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t5-t9+t10-t36)*2.4E+1;
    double t40 = t38*t38;
    double t42 = t41-3.0E+1;
    double t46 = t43+2.0;
    double t47 = t45+1.0;
    double t53 = t29*t38*t48*(3.6E+1/5.0);
    double t50 = (t27*t42)/1.0E+1;
    double t55 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*-2.0;
    double t56 = t27*t42*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*2.0;
    double t57 = t28*t42*t48*(t5-t9+t10-t36)*(-3.0/5.0);
    double t58 = t28*t42*t48*(t5-t9+t10-t36)*(3.0/5.0);
    double t51 = -t50;
    double t59 = t27*t48*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(t5-t9+t10-t36)*6.0;
    double t60 = t28*t38*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*(3.0/5.0);
    double t54 = t44+t51;
    return (q0+Tqd0/5.0)*(t52+t56+t59+pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1)+t2*(t53+t58+t60+t30*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0)*(9.0/2.5E+1)+(t29*((t27*t42)/5.0+t28*(t5-t9+t10-t36)*(1.2E+1/5.0))*pow(t5-t9+t10-t36,3.0))/1.0E+2+t29*t38*t42*((t27*(t5-t9+t10-t36))/5.0-2.0)*(3.0/1.0E+2))-t2*(t31*t40*(3.0/1.25E+2)+t30*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0)*(1.8E+1/2.5E+1)+(t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2+(t29*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*pow(t5-t9+t10-t36,3.0))/5.0E+1+t29*t38*t42*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*(3.0/5.0E+1))-(q0+Tqd0*(2.0/5.0)+TTqdd0/2.0E+1)*(t52+t53+t56+t58+t59+t60)+(t29*((t27*(t5-t9+t10-t36))/5.0-2.0)*pow(t5-t9+t10-t36,3.0))/1.0E+2-(t29*((t27*(t5-t9+t10-t36))/1.0E+1-1.0)*pow(t5-t9+t10-t36,3.0))/5.0E+1+t2*t31*t40*(3.0/1.25E+2)+t28*t38*t48*(3.0/1.0E+1)-q0*pow((t27*(t5-t9+t10-t36))/1.0E+1-1.0,3.0)*(t50+t28*(t5-t9+t10-t36)*(6.0/5.0))*2.0E+1+(t2*t30*t42*pow(t5-t9+t10-t36,3.0))/5.0E+2;
}

double q_des_k_indep(double q0, double Tqd0, double TTqdd0, double s) {
    return q0 + Tqd0*s - 6*Tqd0*pow(s,3) + 8*Tqd0*pow(s,4) - 3*Tqd0*pow(s,5) + (TTqdd0*pow(s,2))*0.5 - (3*TTqdd0*pow(s,3))*0.5 + (3*TTqdd0*pow(s,4))*0.5 - (TTqdd0*pow(s,5))*0.5;
}

double qd_des_k_indep(double q0, double Tqd0, double TTqdd0, double s) {
    return (pow(s - 1,2)*(2*Tqd0 + 4*Tqd0*s + 2*TTqdd0*s - 30*Tqd0*pow(s,2) - 5*TTqdd0*pow(s,2)))*0.5 / DURATION;
}

double qdd_des_k_indep(double q0, double Tqd0, double TTqdd0, double s) {
    return -(s - 1.0)*(TTqdd0 - (36*Tqd0 + 8*TTqdd0)*s + (60*Tqd0 + 10*TTqdd0)*pow(s, 2))  / (DURATION * DURATION);
}

#endif