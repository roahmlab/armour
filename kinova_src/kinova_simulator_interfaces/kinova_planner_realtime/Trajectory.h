#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "PZsparse.h"

// These values are specifically corresponded with a Bezuer curve parameterization
#define QDD_DES_K_DEP_MAXIMA (0.5 - sqrt(3) / 6)
#define QDD_DES_K_DEP_MINIMA (0.5 + sqrt(3) / 6)

// 5th order Bezier curve
// The initial position/velocity/acceleration is equal to q0/qd0/qdd0
// The end position is equal to q0 + k
// The end velocity/acceleration is equal to 0
//
// NOTE:
// This is just a simplified implementation!!!
// t is automatically set to range from 0 to 1, so you don't have to scale.
// Everything has to be changed if the range of t is not [0,1]
//
// -1 <= k <= 1
// k_actual = k * k_range
//
// q_des   = t^3*(6*t^2 - 15*t + 10) * k_actual + 
//           q0 + qd0*t - 6*qd0*t^3 + 8*qd0*t^4 - 3*qd0*t^5 + (qdd0*t^2)/2 - (3*qdd0*t^3)/2 + (3*qdd0*t^4)/2 - (qdd0*t^5)/2
//
// qd_des  = 30*t^2*(t - 1)^2 * k_actual + 
//           ((t - 1)^2*(2*qd0 + 4*qd0*t + 2*qdd0*t - 30*qd0*t^2 - 5*qdd0*t^2))/2
//
// qdd_des = 60*t*(2*t^2 - 3*t + 1) * k_actual + 
//           -(t - 1)*(qdd0 - 36*qd0*t - 8*qdd0*t + 60*qd0*t^2 + 10*qdd0*t^2)
//
class BezierCurve{
public:
    Eigen::VectorXd q0;
    Eigen::VectorXd qd0;
    Eigen::VectorXd qdd0;

    Eigen::VectorXd Tqd0; // qd0 * T
    Eigen::VectorXd TTqdd0; // qdd0 * T ^ 2

    double q_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    double q_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    double q_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    double q_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    double qd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    double qd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    double qd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    double qd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    double qdd_des_k_indep_extrema_1[NUM_FACTORS] = {0.0};
    double qdd_des_k_indep_extrema_2[NUM_FACTORS] = {0.0};
    double qdd_des_k_indep_extremum_1[NUM_FACTORS] = {0.0};
    double qdd_des_k_indep_extremum_2[NUM_FACTORS] = {0.0};

    double ds = 0;

    // rotation matrix (and its transpose) of each joint
    PZsparseArray cos_q_des;
    PZsparseArray sin_q_des;
    PZsparseArray R;
    PZsparseArray R_t;

    // joint velocity
    PZsparseArray qd_des;

    // auxiliary joint velocity
    PZsparseArray qda_des;

    // joint acceleration
    PZsparseArray qdda_des;

    BezierCurve();

    BezierCurve(const Eigen::VectorXd& q0_inp, 
                const Eigen::VectorXd& qd0_inp, 
                const Eigen::VectorXd& qdd0_inp);

    ~BezierCurve() {};

    // convert to polynomial zonotope using 1st/2nd order Taylor expansion
    void makePolyZono(int t_ind);

    // return the min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremum(double* extremum, const double* k) const;

    // return the gradient of min and max of the joint position throughout the whole desired trajectory
    void returnJointPositionExtremumGradient(double* extremumGradient, const double* k) const;

    // return the min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremum(double* extremum, const double* k) const;

    // return the gradient of min and max of the joint velocity throughout the whole desired trajectory
    void returnJointVelocityExtremumGradient(double* extremumGradient, const double* k) const;
};

// helper functions
// q0, qd0, qdd0, k here are scalars since all joints are using the same Bezier curve representation
double q_des_func(double q0, double Tqd0, double TTqdd0, double k, double t);

double qd_des_func(double q0, double Tqd0, double TTqdd0, double k, double t);

double qdd_des_func(double q0, double Tqd0, double TTqdd0, double k, double t);

// derivative of the second extrema of q_des (when qd_des = 0) w.r.t k 
double q_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// derivative of the third extrema of q_des (when qd_des = 0) w.r.t k 
double q_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// derivative of the second extrema of qd_des (when qdd_des = 0) w.r.t k 
double qd_des_extrema2_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// derivative of the third extrema of qd_des (when qdd_des = 0) w.r.t k 
double qd_des_extrema3_k_derivative(double q0, double Tqd0, double TTqdd0, double k);

// k-independent part of q_des
double q_des_k_indep(double q0, double Tqd0, double TTqdd0, double t);

// k-independent part of qd_des
double qd_des_k_indep(double q0, double Tqd0, double TTqdd0, double t);

// k-independent part of qdd_des
double qdd_des_k_indep(double q0, double Tqd0, double TTqdd0, double t);

#endif