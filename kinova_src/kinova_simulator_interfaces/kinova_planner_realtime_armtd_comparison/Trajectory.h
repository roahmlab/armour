#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "PZsparse.h"

// A simple quadratic trajectory whose acceleration is parameterized
//
// -1 <= k <= 1
// k_actual = k * k_range
//
// q_des   = q0 + qd0 * t + 0.5 * k_actual * t^2
//
// qd_des  = qd0 + k_actual * t
//
// qdd_des = k_actual
//
class ConstantAccelerationCurve{
public:
    double* q0 = nullptr;
    double* qd0 = nullptr;

    double dt;

    double* c_cos_q_des = nullptr; // center of zonotope
    double* g_cos_q_des = nullptr; // k-dependent generator of zonotope
    double* r_cos_q_des = nullptr; // radius of zonotope

    double* c_sin_q_des = nullptr; // center of zonotope
    double* g_sin_q_des = nullptr; // k-dependent generator of zonotope
    double* r_sin_q_des = nullptr; // radius of zonotope

    // PZsparse cos_q_des[NUM_TIME_STEPS * NUM_FACTORS];
    // PZsparse sin_q_des[NUM_TIME_STEPS * NUM_FACTORS];

    PZsparseArray R;
    PZsparseArray R_t;

    double* k_range = nullptr;

    ConstantAccelerationCurve() {};

    ConstantAccelerationCurve(double* q0_inp, double* qd0_inp, 
                              double* c_cos_q_des_inp, double* g_cos_q_des_inp, double* r_cos_q_des_inp,
                              double* c_sin_q_des_inp, double* g_sin_q_des_inp, double* r_sin_q_des_inp,
                              double* k_range_inp);

    ~ConstantAccelerationCurve() {};

    // convert to polynomial zonotope using 1st/2nd order Taylor expansion
    void makePolyZono(int t_ind);

    // return the min and max of the joint position and velocity throughout the whole desired trajectory
    void returnJointStateExtremum(double* extremum, const double* k) const;

    // return the graident of the min and max of the joint position and velocity throughout the whole desired trajectory
    void returnJointStateExtremumGradient(double* extremumGradient, const double* k) const;
};

#endif