#ifndef ROBOT_MODEL_HPP
#define ROBOT_MODEL_HPP

#include "spatial.cpp"
#include "spatial_interval.cpp"

class Model {
public:
    // Spanning tree structure variables
    int numJoints; // Number of joints in tree
    std::vector<Twist> S; // Joint twists
    std::vector<RigidInertia> I; // Link spatial inertias
    std::vector<Transform> XTree; //Body tree transforms
    std::vector<int> lam; // Parent array that details joint connectivity
    std::vector<Transform> CoM; // CoM transforms
    Eigen::VectorXd transI; // Transmission inertias
    Eigen::VectorXd gearRatios;
    Twist gravity;

    // Additional elements from MuJoCo
    Eigen::VectorXd friction;
    Eigen::VectorXd damping;
    
    Model(int numJoints_in = 20);

    Model(std::string filename);
};

class IntModel {
public:
    // Spanning tree structure variables
    int numJoints; // Number of joints in tree
    std::vector<IntTwist> S; // Joint twists
    std::vector<IntRigidInertia> I; // Link spatial inertias
    std::vector<IntTransform> XTree; // Body tree transforms
    std::vector<int> lam; // Parent array that details joint connectivity
    std::vector<IntTransform> CoM; // CoM transforms
    Eigen::Matrix<Interval, Eigen::Dynamic, 1> transI; // Transmission inertias
    Eigen::VectorXd gearRatios;
    IntTwist gravity;

    // Additional elements from MuJoCo
    Eigen::VectorXd friction;
    Eigen::VectorXd damping;

    // Storage variables for Interval RNEA
    // Declare velocity, acceleration and force arrays
    std::vector<IntTwist> v;
    std::vector<IntTwist> va;
    std::vector<IntTwist> a;
    std::vector<IntWrench> f;
    std::vector<IntTransform> Xbw; //Transform from body to world frame
    std::vector<IntTwist> Sb; //Screw axis in body frame
    std::vector<IntTransform> Xli_i; //Transform from body i-1 to i

    IntModel(int numJoints_in = 20);

    IntModel(const Model* model, double eps);
};

class Robot {
public:
    int numJoints = 0;
    Model* RobotModelPtr = nullptr;
    IntModel* IntRobotModelPtr = nullptr;

    bool ifConstrained = false;

    explicit Robot() {};

    explicit Robot(std::string filename, double eps);

    ~Robot();

    // Define the constraint function if the robot is a constrained system.
    // Fill in the constrained torque for unactuated dimensions in tau.
    virtual void applyConstraints(const Eigen::VectorXd &q, Eigen::VectorXd &tau) {}
    virtual void applyConstraints_Int(const Eigen::VectorXd &q, VectorXint &tau) {}

    // Compute integration error.
    // Could be different if the robot is constrained.
    virtual double integrateStateError(const Eigen::VectorXd &q, const Eigen::VectorXd &qd);

    // An virtual function for extra procedures to compute torque
    virtual Eigen::VectorXd finalizeTorque(Eigen::VectorXd& tau);

    // Virtual function for filling in unactuated joints
    virtual void fillUnactJoints(Eigen::VectorXd &q, Eigen::VectorXd &q_d, Eigen::VectorXd &q_dd);
};

#endif