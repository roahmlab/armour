#ifndef TRAJECTORY_CPP
#define TRAJECTORY_CPP

#include "Trajectory.h"

ConstantAccelerationCurve::ConstantAccelerationCurve(double* q0_inp, double* qd0_inp, 
                                                     double* c_cos_q_des_inp, double* g_cos_q_des_inp, double* r_cos_q_des_inp,
                                                     double* c_sin_q_des_inp, double* g_sin_q_des_inp, double* r_sin_q_des_inp,
                                                     double* k_range_inp) {
    q0 = q0_inp;
    qd0 = qd0_inp; 

    c_cos_q_des = c_cos_q_des_inp;
    g_cos_q_des = g_cos_q_des_inp;
    r_cos_q_des = r_cos_q_des_inp;

    c_sin_q_des = c_sin_q_des_inp;
    g_sin_q_des = g_sin_q_des_inp;
    r_sin_q_des = r_sin_q_des_inp;

    R = PZsparseArray(NUM_JOINTS + 1, NUM_TIME_STEPS);
    R_t = PZsparseArray(NUM_JOINTS, NUM_TIME_STEPS);

    k_range = k_range_inp;

    dt = 1.0 / NUM_TIME_STEPS;
}

void ConstantAccelerationCurve::makePolyZono(int t_ind) {
    assert(t_ind < NUM_TIME_STEPS);

    for (int i = 0; i < NUM_FACTORS; i++) {
        const double k_range_elt = k_range[i];
        const double cos_q0 = cos(q0[i]);
        const double sin_q0 = sin(q0[i]);

        // cos(q_des)
        double cos_q_des_center = cos_q0 * c_cos_q_des[i * NUM_TIME_STEPS + t_ind] - sin_q0 * c_sin_q_des[i * NUM_TIME_STEPS + t_ind];
        double cos_q_des_coeff[2];
        cos_q_des_coeff[0] = cos_q0 * g_cos_q_des[i * NUM_TIME_STEPS + t_ind] - sin_q0 * g_sin_q_des[i * NUM_TIME_STEPS + t_ind];
        cos_q_des_coeff[1] = fabs(cos_q0) * r_cos_q_des[i * NUM_TIME_STEPS + t_ind] + fabs(sin_q0) * r_sin_q_des[i * NUM_TIME_STEPS + t_ind];

        cos_q_des_coeff[1] *= 4.0;

        uint64_t cos_q_des_degree[2][NUM_FACTORS * 6] = {0};
        cos_q_des_degree[0][i] = 1; // k
        cos_q_des_degree[1][i + NUM_FACTORS * 4] = 1; // cosqe
        // cos_q_des[t_ind * NUM_FACTORS + i] = PZsparse(cos_q_des_center, cos_q_des_coeff, cos_q_des_degree, 2);

        // sin(q_des)
        double sin_q_des_center = cos_q0 * c_sin_q_des[i * NUM_TIME_STEPS + t_ind] + sin_q0 * c_cos_q_des[i * NUM_TIME_STEPS + t_ind];
        double sin_q_des_coeff[2];
        sin_q_des_coeff[0] = cos_q0 * g_sin_q_des[i * NUM_TIME_STEPS + t_ind] + sin_q0 * g_cos_q_des[i * NUM_TIME_STEPS + t_ind];
        sin_q_des_coeff[1] = fabs(cos_q0) * r_sin_q_des[i * NUM_TIME_STEPS + t_ind] + fabs(sin_q0) * r_cos_q_des[i * NUM_TIME_STEPS + t_ind];

        sin_q_des_coeff[1] *= 4.0;

        uint64_t sin_q_des_degree[2][NUM_FACTORS * 6] = {0};
        sin_q_des_degree[0][i] = 1; // k
        sin_q_des_degree[1][i + NUM_FACTORS * 5] = 1; // sinqe
        // sin_q_des[t_ind * NUM_FACTORS + i] = PZsparse(sin_q_des_center, sin_q_des_coeff, sin_q_des_degree, 2);

        R(i, t_ind) = PZsparse(rots[i * 3], rots[i * 3 + 1], rots[i * 3 + 2]);

        if (axes[i] != 0) {
            R(i, t_ind) = R(i, t_ind) * PZsparse(cos_q_des_center, cos_q_des_coeff, cos_q_des_degree, 2,
                                                 sin_q_des_center, sin_q_des_coeff, sin_q_des_degree, 2,
                                                 axes[i]);
        }

        R_t(i, t_ind) = R(i, t_ind).transpose();
    }

    // assume all fixed joints are at the end of the kinematics chain
    for (int i = NUM_FACTORS; i < NUM_JOINTS; i++) {
        R(i, t_ind) = PZsparse(rots[i * 3], rots[i * 3 + 1], rots[i * 3 + 2]);
        R_t(i, t_ind) = R(i, t_ind).transpose();
    }

    R(NUM_JOINTS, t_ind) = PZsparse(0, 0, 0);
}

void ConstantAccelerationCurve::returnJointStateExtremum(double* extremum, const double* k) const {
    double t_move = 0.5;
    double t_total = 1.0;
    double t_to_stop = t_total - t_move;

    for (int i = 0; i < NUM_FACTORS; i++){
        double k_actual = k_range[i] * k[i];
        double q_peak = q0[i] + qd0[i] * t_move + k_actual * t_move * t_move * 0.5;
        double q_dot_peak = qd0[i] + k_actual * t_move;
        double q_ddot_to_stop = -q_dot_peak / t_to_stop;
        double q_stop = q_peak + q_dot_peak * t_to_stop + 0.5 * q_ddot_to_stop * t_to_stop * t_to_stop;
        double t_max_min_to_peak = -qd0[i] / k_actual;

        double q_max_to_peak;
        double q_min_to_peak;
        double q_dot_max_to_peak;
        double q_dot_min_to_peak;

        double grad_q_max_to_peak;
        double grad_q_min_to_peak;
        double grad_q_dot_max_to_peak;
        double grad_q_dot_min_to_peak;

        double q_max_to_stop;
        double q_min_to_stop;
        double q_dot_max_to_stop;
        double q_dot_min_to_stop;

        double grad_q_max_to_stop;
        double grad_q_min_to_stop;
        double grad_q_dot_max_to_stop;
        double grad_q_dot_min_to_stop;

        double q_endpoints_ordered[2];
        double grad_q_endpoints_ordered[2];

        if (q_peak >= q0[i]){
            q_endpoints_ordered[0] = q0[i]; 
            q_endpoints_ordered[1] = q_peak;
            grad_q_endpoints_ordered[0] = 0; 
            grad_q_endpoints_ordered[1] = 0.5 * t_move * t_move;
        }
        else{
            q_endpoints_ordered[0] = q_peak; 
            q_endpoints_ordered[1] = q0[i];
            grad_q_endpoints_ordered[0] = 0.5 * t_move * t_move; 
            grad_q_endpoints_ordered[1] = 0;
        }

        if (t_max_min_to_peak > 0 && t_max_min_to_peak < t_move){
            if (k_actual >= 0){
                q_min_to_peak = q0[i] + qd0[i] * t_max_min_to_peak + 0.5 * k_actual * t_max_min_to_peak * t_max_min_to_peak;
                q_max_to_peak  = q_endpoints_ordered[1];
                grad_q_min_to_peak = (0.5 * qd0[i] * qd0[i]) / (k_actual * k_actual);
                grad_q_max_to_peak = grad_q_endpoints_ordered[1];
            }
            else{
                q_min_to_peak = q_endpoints_ordered[0];
                q_max_to_peak = q0[i] + qd0[i] * t_max_min_to_peak + 0.5 * k_actual * t_max_min_to_peak * t_max_min_to_peak;
                grad_q_min_to_peak = grad_q_endpoints_ordered[0];
                grad_q_max_to_peak = (0.5 * qd0[i] * qd0[i]) / (k_actual * k_actual);
            }
        }
        else{
            q_min_to_peak = q_endpoints_ordered[0];
            q_max_to_peak = q_endpoints_ordered[1];
            
            grad_q_min_to_peak = grad_q_endpoints_ordered[0];
            grad_q_max_to_peak = grad_q_endpoints_ordered[1];
        }

        if( q_dot_peak >= qd0[i]){
            q_dot_min_to_peak = qd0[i];
            q_dot_max_to_peak = q_dot_peak;
            
            grad_q_dot_min_to_peak = 0;
            grad_q_dot_max_to_peak = t_move;
        }
        else{
            q_dot_min_to_peak = q_dot_peak;
            q_dot_max_to_peak = qd0[i];
            
            grad_q_dot_min_to_peak = t_move;
            grad_q_dot_max_to_peak = 0;
        }

        if( q_stop >= q_peak){
            q_min_to_stop = q_peak;
            q_max_to_stop = q_stop;
            
            grad_q_min_to_stop = 0.5 * t_move * t_move;
            grad_q_max_to_stop = 0.5 * t_move * t_move + 0.5 * t_move * t_to_stop;
        }
        else{
            q_min_to_stop = q_stop;
            q_max_to_stop = q_peak;
            
            grad_q_min_to_stop = 0.5 * t_move * t_move + 0.5 * t_move * t_to_stop;
            grad_q_max_to_stop = 0.5 * t_move * t_move;
        }

        if(q_dot_peak >= 0){
            q_dot_min_to_stop = 0;
            q_dot_max_to_stop = q_dot_peak;
            
            grad_q_dot_min_to_stop = 0;
            grad_q_dot_max_to_stop = t_move;
        }
        else{
            q_dot_min_to_stop = q_dot_peak;
            q_dot_max_to_stop = 0;
            
            grad_q_dot_min_to_stop = t_move;
            grad_q_dot_max_to_stop = 0;
        }

        if (q_min_to_peak <= q_min_to_stop){
            extremum[i] = q_min_to_peak; // q_min[i]
        }
        else{
            extremum[i] = q_min_to_stop; // q_min[i]
        }

        if (q_max_to_peak >= q_max_to_stop){
            extremum[i + NUM_FACTORS] = q_max_to_peak; // q_max[i]
        }
        else{
            extremum[i + NUM_FACTORS] = q_max_to_stop; // q_max[i]
        }

        if (q_dot_min_to_peak <= q_dot_min_to_stop){
            extremum[i + 2 * NUM_FACTORS] = q_dot_min_to_peak; // q_dot_min[i]
        }
        else{
            extremum[i + 2 * NUM_FACTORS] = q_dot_min_to_stop; // q_dot_min[i]
        }

        if (q_dot_max_to_peak >= q_dot_max_to_stop){
            extremum[i + 3 * NUM_FACTORS] = q_dot_max_to_peak; // q_dot_max[i]
        }
        else{
            extremum[i + 3 * NUM_FACTORS] = q_dot_max_to_stop; // q_dot_max[i]
        }
    }
}

void ConstantAccelerationCurve::returnJointStateExtremumGradient(double* extremumGradient, const double* k) const {
    memset(extremumGradient, 0, 4 * NUM_FACTORS * NUM_FACTORS);

    double t_move = 0.5;
    double t_total = 1.0;
    double t_to_stop = t_total - t_move;

    for (int i = 0; i < NUM_FACTORS; i++){
        double k_actual = k_range[i] * k[i];
        double q_peak = q0[i] + qd0[i] * t_move + k_actual * t_move * t_move * 0.5;
        double q_dot_peak = qd0[i] + k_actual * t_move;
        double q_ddot_to_stop = -q_dot_peak / t_to_stop;
        double q_stop = q_peak + q_dot_peak * t_to_stop + 0.5 * q_ddot_to_stop * t_to_stop * t_to_stop;
        double t_max_min_to_peak = -qd0[i] / k_actual;

        double q_max_to_peak;
        double q_min_to_peak;
        double q_dot_max_to_peak;
        double q_dot_min_to_peak;

        double grad_q_max_to_peak;
        double grad_q_min_to_peak;
        double grad_q_dot_max_to_peak;
        double grad_q_dot_min_to_peak;

        double q_max_to_stop;
        double q_min_to_stop;
        double q_dot_max_to_stop;
        double q_dot_min_to_stop;

        double grad_q_max_to_stop;
        double grad_q_min_to_stop;
        double grad_q_dot_max_to_stop;
        double grad_q_dot_min_to_stop;

        double q_endpoints_ordered[2];
        double grad_q_endpoints_ordered[2];

        if (q_peak >= q0[i]){
            q_endpoints_ordered[0] = q0[i]; 
            q_endpoints_ordered[1] = q_peak;
            grad_q_endpoints_ordered[0] = 0; 
            grad_q_endpoints_ordered[1] = 0.5 * t_move * t_move;
        }
        else{
            q_endpoints_ordered[0] = q_peak; 
            q_endpoints_ordered[1] = q0[i];
            grad_q_endpoints_ordered[0] = 0.5 * t_move * t_move; 
            grad_q_endpoints_ordered[1] = 0;
        }

        if (t_max_min_to_peak > 0 && t_max_min_to_peak < t_move){
            if (k_actual >= 0){
                q_min_to_peak = q0[i] + qd0[i] * t_max_min_to_peak + 0.5 * k_actual * t_max_min_to_peak * t_max_min_to_peak;
                q_max_to_peak  = q_endpoints_ordered[1];
                grad_q_min_to_peak = (0.5 * qd0[i] * qd0[i]) / (k_actual * k_actual);
                grad_q_max_to_peak = grad_q_endpoints_ordered[1];
            }
            else{
                q_min_to_peak = q_endpoints_ordered[0];
                q_max_to_peak = q0[i] + qd0[i] * t_max_min_to_peak + 0.5 * k_actual * t_max_min_to_peak * t_max_min_to_peak;
                grad_q_min_to_peak = grad_q_endpoints_ordered[0];
                grad_q_max_to_peak = (0.5 * qd0[i] * qd0[i]) / (k_actual * k_actual);
            }
        }
        else{
            q_min_to_peak = q_endpoints_ordered[0];
            q_max_to_peak = q_endpoints_ordered[1];
            
            grad_q_min_to_peak = grad_q_endpoints_ordered[0];
            grad_q_max_to_peak = grad_q_endpoints_ordered[1];
        }

        if( q_dot_peak >= qd0[i]){
            q_dot_min_to_peak = qd0[i];
            q_dot_max_to_peak = q_dot_peak;
            
            grad_q_dot_min_to_peak = 0;
            grad_q_dot_max_to_peak = t_move;
        }
        else{
            q_dot_min_to_peak = q_dot_peak;
            q_dot_max_to_peak = qd0[i];
            
            grad_q_dot_min_to_peak = t_move;
            grad_q_dot_max_to_peak = 0;
        }

        if( q_stop >= q_peak){
            q_min_to_stop = q_peak;
            q_max_to_stop = q_stop;
            
            grad_q_min_to_stop = 0.5 * t_move * t_move;
            grad_q_max_to_stop = 0.5 * t_move * t_move + 0.5 * t_move * t_to_stop;
        }
        else{
            q_min_to_stop = q_stop;
            q_max_to_stop = q_peak;
            
            grad_q_min_to_stop = 0.5 * t_move * t_move + 0.5 * t_move * t_to_stop;
            grad_q_max_to_stop = 0.5 * t_move * t_move;
        }

        if(q_dot_peak >= 0){
            q_dot_min_to_stop = 0;
            q_dot_max_to_stop = q_dot_peak;
            
            grad_q_dot_min_to_stop = 0;
            grad_q_dot_max_to_stop = t_move;
        }
        else{
            q_dot_min_to_stop = q_dot_peak;
            q_dot_max_to_stop = 0;
            
            grad_q_dot_min_to_stop = t_move;
            grad_q_dot_max_to_stop = 0;
        }

        if (q_min_to_peak <= q_min_to_stop){
            // extremum[i] = q_min_to_peak; // q_min[i]
            extremumGradient[i * NUM_FACTORS + i] = grad_q_min_to_peak;
        }
        else{
            // extremum[i] = q_min_to_stop; // q_min[i]
            extremumGradient[i * NUM_FACTORS + i] = grad_q_min_to_stop;
        }

        if (q_max_to_peak >= q_max_to_stop){
            // extremum[i + NUM_FACTORS] = q_max_to_peak; // q_max[i]
            extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + i] = grad_q_max_to_peak;
        }
        else{
            // extremum[i + NUM_FACTORS] = q_max_to_stop; // q_max[i]
            extremumGradient[(i + NUM_FACTORS) * NUM_FACTORS + i] = grad_q_max_to_stop;
        }

        if (q_dot_min_to_peak <= q_dot_min_to_stop){
            // extremum[i + 2 * NUM_FACTORS] = q_dot_min_to_peak; // q_dot_min[i]
            extremumGradient[(i + 2 * NUM_FACTORS) * NUM_FACTORS + i] = grad_q_dot_min_to_peak;
        }
        else{
            // extremum[i + 2 * NUM_FACTORS] = q_dot_min_to_stop; // q_dot_min[i]
            extremumGradient[(i + 2 * NUM_FACTORS) * NUM_FACTORS + i] = grad_q_dot_min_to_stop;
        }

        if (q_dot_max_to_peak >= q_dot_max_to_stop){
            // extremum[i + 3 * NUM_FACTORS] = q_dot_max_to_peak; // q_dot_max[i]
            extremumGradient[(i + 3 * NUM_FACTORS) * NUM_FACTORS + i] = grad_q_dot_max_to_peak;
        }
        else{
            // extremum[i + 3 * NUM_FACTORS] = q_dot_max_to_stop; // q_dot_max[i]
            extremumGradient[(i + 3 * NUM_FACTORS) * NUM_FACTORS + i] = grad_q_dot_max_to_stop;
        }
    }
}

#endif