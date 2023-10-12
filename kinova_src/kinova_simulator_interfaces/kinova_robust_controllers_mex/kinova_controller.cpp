#include <memory>
#include "mex.h"
#include "robust_controller.hpp"
#include "robot_models.hpp"
#include "robot_model_file_path.hpp"

void copyPointerToEigenVector(Eigen::VectorXd& a, double* b, unsigned int n) {
    for (unsigned int i = 0; i < n; i++) {
        a(i) = b[i];
    }
}

void copyEigenVectorToPointer(double* a, Eigen::VectorXd& b, unsigned int n) {
    for (unsigned int i = 0; i < n; i++) {
        a[i] = b(i);
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // create a Kinova model
    double ModelUncertainty = 0.03;
    if (nrhs > 9) {
        ModelUncertainty = *(double*)mxGetData(prhs[9]);
    }

    std::unique_ptr<Robot> Kinova = nullptr;
    try {
        Kinova = std::make_unique<Robot>(RobotFilePath, ModelUncertainty);
    }
    catch (int error_code) {
        mexErrMsgTxt("Wrong model file path! Check armour-dev/kinova_src/kinova_simulator_interfaces/kinova_robust_controllers_mex/kinova_controller.cpp!");
    }

    int numJoints = Kinova->numJoints;

    double* Kr_input = (double*)mxGetData(prhs[0]);
    Eigen::MatrixXd Kr = Eigen::MatrixXd::Identity(numJoints, numJoints);
    for (int i = 0; i < numJoints; i++) {
        Kr(i,i) = Kr_input[i];
    }

    double* alpha_input = (double*)mxGetData(prhs[1]);

    double* V_max_input = (double*)mxGetData(prhs[2]);

    double* r_norm_threshold_input = (double*)mxGetData(prhs[3]);
    
    // create a robust controller
    RobustController c(Kr, *alpha_input, *V_max_input, *r_norm_threshold_input);
    c.applyFriction = false;

    double* q_input = (double*)mxGetData(prhs[4]);
    Eigen::VectorXd q(numJoints);
    copyPointerToEigenVector(q, q_input, numJoints);

    double* q_d_input = (double*)mxGetData(prhs[5]);
    Eigen::VectorXd q_d(numJoints);
    copyPointerToEigenVector(q_d, q_d_input, numJoints);

    double* qd_input = (double*)mxGetData(prhs[6]);
    Eigen::VectorXd qd(numJoints);
    copyPointerToEigenVector(qd, qd_input, numJoints);

    double* qd_d_input = (double*)mxGetData(prhs[7]);
    Eigen::VectorXd qd_d(numJoints);
    copyPointerToEigenVector(qd_d, qd_d_input, numJoints);

    double* qd_dd_input = (double*)mxGetData(prhs[8]);
    Eigen::VectorXd qd_dd(numJoints);
    copyPointerToEigenVector(qd_dd, qd_dd_input, numJoints);

    Eigen::VectorXd u = c.update(Kinova.get(), q, q_d, qd, qd_d, qd_dd);

    plhs[0] = mxCreateNumericMatrix(numJoints, 1, mxDOUBLE_CLASS, mxREAL);
    copyEigenVectorToPointer((double*)mxGetData(plhs[0]), u, numJoints);

    plhs[1] = mxCreateNumericMatrix(numJoints, 1, mxDOUBLE_CLASS, mxREAL);
    copyEigenVectorToPointer((double*)mxGetData(plhs[1]), c.u_nominal, numJoints);

    plhs[2] = mxCreateNumericMatrix(numJoints, 1, mxDOUBLE_CLASS, mxREAL);
    copyEigenVectorToPointer((double*)mxGetData(plhs[2]), c.v, numJoints);

    nlhs = 3;
}
