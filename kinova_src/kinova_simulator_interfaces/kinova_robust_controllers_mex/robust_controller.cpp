#ifndef ROBUST_CONTROLLER_CPP
#define ROBUST_CONTROLLER_CPP

#include "robust_controller.hpp"

RobustController::RobustController() {
    Kr = Eigen::VectorXd::Zero(1);
    Kp = Eigen::VectorXd::Zero(2);
    Ki = Eigen::VectorXd::Zero(2);
    maxErrorBound = 0;
    ssErrorAccum = 0;

    alpha = 0;
    V_max = 0;
    q_d_zero = Eigen::VectorXd::Zero(1);
    qa_d_zero = Eigen::VectorXd::Zero(1);

    u_nominal = Eigen::VectorXd::Zero(1);
    v = Eigen::VectorXd::Zero(1);
}

RobustController::RobustController(Eigen::MatrixXd Kr_in, Eigen::Vector2d Kp_in, Eigen::Vector2d Ki_in, double maxError_in) {
    // if (which_method != ROBUST_INPUT_METHOD::ALTHOFF) {
    //     cout << "Wrong constructor for this robust input method!\n";
    //     throw;
    // }
    which_method = ROBUST_INPUT_METHOD::ALTHOFF;

    Kr = Kr_in;
    Kp = Kp_in;
    Ki = Ki_in;
    maxErrorBound = maxError_in;

    // Initialize state variables
    ssErrorAccum = 0;

    int numJoints = Kr.rows();
    u_nominal = Eigen::VectorXd::Zero(numJoints);
    v = Eigen::VectorXd::Zero(numJoints);
}

RobustController::RobustController(Eigen::MatrixXd Kr_in, double alpha_in, double V_max_in, double r_norm_threshold_in) {
    // if (which_method != ROBUST_INPUT_METHOD::ARMOUR) {
    //     cout << "Wrong constructor for this robust input method!\n";
    //     throw;
    // }
    which_method = ROBUST_INPUT_METHOD::ARMOUR;

    Kr = Kr_in;
    alpha = alpha_in;
    V_max = V_max_in;
    r_norm_threshold = r_norm_threshold_in;

    // Initialize zero vectors
    int numJoints = Kr.rows();
    q_d_zero = Eigen::VectorXd::Zero(numJoints);
    qa_d_zero = Eigen::VectorXd::Zero(numJoints);

    u_nominal = Eigen::VectorXd::Zero(numJoints);
    v = Eigen::VectorXd::Zero(numJoints);
}

Eigen::VectorXd RobustController::update(Robot* RobotPtr,
                        Eigen::VectorXd &q,  Eigen::VectorXd &q_d, 
                        Eigen::VectorXd &qd, Eigen::VectorXd &qd_d, Eigen::VectorXd &qd_dd, 
                        double deltaT, double eAcc) {
    int numJoints = RobotPtr->numJoints;

    // Step 1, calculate reference terms
    Eigen::VectorXd q_diff = qd - q;

    // clamp joint positions to [-pi,pi) for continuous revolute joints
    for (int i = 0; i < numJoints; i++) {
        q_diff(i) = clamp(q_diff(i));
    }
    
    Eigen::VectorXd qa_d = qd_d + Kr*q_diff;
    Eigen::VectorXd qa_dd = qd_dd + Kr*(qd_d - q_d);
    Eigen::VectorXd r = qd_d - q_d + Kr*q_diff;

    // Step 2, calculate nominal torque from passivity RNEA
    passRNEA(u_nominal, RobotPtr->RobotModelPtr, q, q_d, qa_d, qa_dd, applyFriction);
    
    // Step 3, calculate interval torque from passivity RNEA
    VectorXint u_interval(numJoints);
    passRNEA_Int(u_interval, RobotPtr->IntRobotModelPtr, q, q_d, qa_d, qa_dd, applyFriction);

    // Step 3.5, apply constraints to get actuated torques if the robot is a constrained model
    if (RobotPtr->ifConstrained) {
        RobotPtr->applyConstraints(q, u_nominal);
        RobotPtr->applyConstraints_Int(q, u_interval);
    }
    
    // Interval check
    for (size_t i = 0; i < numJoints; i++) {
        if (u_nominal[i] > u_interval[i].upper() || u_nominal[i] < u_interval[i].lower()) {
            cout << "robust_controller.cpp: update(): Nominal model output falls outside interval output!!!" << endl;
            cout << u_nominal[i] << " " << u_interval[i].lower() << " " << u_interval[i].upper() << endl;
            throw;
        }
    }

    // Step 4, calculate error between nominal and interval
    // Because of constraints, unactuated terms should be 0
    VectorXint phi = u_interval - u_nominal.cast<Interval>();

    // Step 5, calculate error bound (only includes actuated terms)
    Eigen::VectorXd bound = calculateErrorBound(phi);

    v = Eigen::VectorXd::Zero(numJoints);

    if (which_method == ROBUST_INPUT_METHOD::ALTHOFF) {
        // Step 6a, calculate y(t), k(t), only includes actuated terms for state error
        if (calcError) {
            // Only do state error for actuated joints
            double stateError = RobotPtr->integrateStateError(q, qd);
            if (stateError > maxErrorBound) {
                ssErrorAccum += stateError * deltaT;
            }
            eAcc = ssErrorAccum;
        }
        
        double phi_t   = Kp[0] + Ki[0]*eAcc;
        double kappa_t = Kp[1] + Ki[1]*eAcc;

        // Step 6b, compute robust input
        v = -(kappa_t*bound.norm() + phi_t) * r;
    }
    else {
        if (RobotPtr->ifConstrained) {
            cout << "ARMOUR robust input method currently does not support constrained models!\n";
            throw;
        }

        double r_norm = r.norm();
        
        // Only compute robust input when r is large enough to avoid numerical issues
        if (r_norm > r_norm_threshold) {
            // Step 6a, calculate V = 0.5 * r' * M * r
            VectorXint Mr(numJoints);
            passRNEA_Int(Mr, RobotPtr->IntRobotModelPtr, q, q_d_zero, qa_d_zero, r, false, false);

            // dot product between r and M * r 
            Interval V_int = 0;
            for (int i = 0; i < numJoints; i++) {
                V_int += 0.5 * r(i) * Mr(i);
            }

            // Step 6b, compute robust input
            double V_sup = V_int.upper();

            double h = -V_sup + V_max;

            double lambda = max(0.0, -alpha * h / r_norm + bound.norm());

            v = -lambda * r / r_norm;
        }
        // else {
        //     v is just zero vector then
        // }
    }

    // Step 7a, compute torque
    Eigen::VectorXd tau = u_nominal - v;

    // Step 7b, customized finalization procedure
    return RobotPtr->finalizeTorque(tau);
}

Eigen::VectorXd RobustController::calculateErrorBound(VectorXint phi) {
    Eigen::VectorXd bound(phi.rows());

    for (int i = 0; i < phi.rows(); i++) {
        bound[i] = max(abs(phi[i].lower()), abs(phi[i].upper()));
    }

    return bound;
}

void RobustController::clearIntegrateError() {
    ssErrorAccum = 0;
}

#endif