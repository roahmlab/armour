#ifndef ROBUST_CONTROLLER_HPP
#define ROBUST_CONTROLLER_HPP

#include "rnea.hpp"
#include "robot_models.hpp"

#define M_TWOPI 6.283185307179586476925286766559

enum class ROBUST_INPUT_METHOD {ALTHOFF, ARMOUR};

inline double clamp(double inp) {
    double res = inp;
    while (res >= M_PI) res -= M_TWOPI;
    while (res < -M_PI) res += M_TWOPI;
    return res;
}

class RobustController {
public:
    // Parameters
    Eigen::MatrixXd Kr;

    ROBUST_INPUT_METHOD which_method = ROBUST_INPUT_METHOD::ALTHOFF;
    //     Method (1): Ultimate Robust Performance Control of Rigid Robot Manipulators using Interval Arithmetic, Andrea Giusti and Matthias Althoff
    Eigen::Vector2d Kp;
    Eigen::Vector2d Ki;
    double maxErrorBound;

    // State variables
    double ssErrorAccum;

    // Can choose not calculate integration error if debug in ode45 in Matlab
    bool calcError = true;

    //     Method (2): ARMOUR
    double alpha = 1.0;
    double V_max = 1e-5;
    double r_norm_threshold = 1e-7;  

    // They are just zero vector place holders for computing M*r
    Eigen::VectorXd q_d_zero;
    Eigen::VectorXd qa_d_zero;

    // Our matlab model does not include friction, so if test with Matlab, turn this off
    bool applyFriction = true;

    // For debug
    Eigen::VectorXd u_nominal;
    Eigen::VectorXd v;
    
    RobustController();

    RobustController(Eigen::MatrixXd Kr_in, Eigen::Vector2d Kp_in, Eigen::Vector2d Ki_in, double maxError_in);

    RobustController(Eigen::MatrixXd Kr_in, double alpha_in, double V_max_in, double r_norm_threshold_in);

    /* Robust Control
    * 
    * Inputs
    *   q - Nx1 joint angle state vector
    *   q_d - Nx1 joint velocity state vector
    *   qd - Nx1 reference joint angle
    *   qd_d - Nx1 reference joint velocity
    *   qd_dd - Nx1 reference joint acceleration
    *
    * Outputs
    *   u - Nx1 joint torque vector
    *   
    * Steps
    * 1. Calculate reference terms qa_d, qa_dd, r
    * 2. u_nominal = passRNEA(q, q_d, qa_d, qa_dd, robot.nominal_model)
    *    This gives the nominal torques required to follow the trajectory
    * 3. u_interval = passRNEA(q, q_d, qa_d, qa_dd, robot.interval_model)
    *    This gives the estimated bound/interval on the torques
    * 4. phi = u_interval - u_nominal - d (leaving d out for now, it is for disturbances)
    *    This is an interval on the error between the nominal and interval torques
    * 5. errorBound = norm(max( abs(phi.lower), abs(phi.upper) ))
    *    Maximum error bound norm
    * ROBUST_INPUT_METHOD::ALTHOFF:
    *   6a. Calculate y(t) and k(t), which are essentially Kp + Ki*integral(state space error over time)
    *   6b. v = - (k(t)*errorBound + y(t)) * (eq_d + Kr*eq)
    * ROBUST_INPUT_METHOD::ARMOUR:
    *   6a. Calculate V = 0.5 * r' * M * r
    *   6b. v = max(0, -alpha * (V_max - V_sup) / ||r|| + errorBound) * r / ||r||
    * 7. u = u_nominal - v
    */
    Eigen::VectorXd update(Robot* RobotPtr,
                           Eigen::VectorXd &q,  Eigen::VectorXd &q_d, 
                           Eigen::VectorXd &qd, Eigen::VectorXd &qd_d, Eigen::VectorXd &qd_dd, 
                           double deltaT = 0, double eAcc = 0);

    Eigen::VectorXd calculateErrorBound(VectorXint phi);

    void clearIntegrateError();
};

#endif