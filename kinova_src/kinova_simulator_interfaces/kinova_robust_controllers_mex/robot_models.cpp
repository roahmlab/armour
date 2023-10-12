#ifndef ROBOT_MODEL_CPP
#define ROBOT_MODEL_CPP

#include "robot_models.hpp"

Model::Model(int numJoints_in) {
    numJoints = numJoints_in;

    S.resize(numJoints);
    I.resize(numJoints);
    XTree.resize(numJoints);
    CoM.resize(numJoints);
    lam.resize(numJoints);
    transI.resize(numJoints);
    gearRatios.resize(numJoints);
    friction.resize(numJoints);
    damping.resize(numJoints);
}

Model::Model(std::string filename) {
    std::ifstream robotFile(filename);
    if (!robotFile.is_open() || filename == "empty") {
        std::cout << "Could not open the robot model file. Check the file path!\n";
        throw 1;
    } 
    else if (robotFile.is_open()) {
        std::string currentRowString;

        while (getline(robotFile, currentRowString)) {
            // temporary variables in string
            std::string field = "";
            std::string tmpInd = ""; 
            std::string tmpValues = "";

            // variables in int and doubles
            int ind;
            std::vector<double> values;

            // helper variables 
            bool notArray = true;
            int i = 0;

            while (i < currentRowString.size()) {
                if (std::isalpha(currentRowString[i]) || currentRowString[i] == '_') {
                    field += currentRowString[i];
                } else if (std::isdigit(currentRowString[i]) && notArray) {
                    tmpInd += currentRowString[i];
                } else if (currentRowString[i] == '<') {
                    notArray = false;
                } else if (currentRowString[i] == '>') {
                    notArray = true;
                    break;
                } else if (!notArray) {
                    tmpValues += currentRowString[i];
                }
                i++;
            }

            // find field and value
            if (tmpInd != "") {
                ind = std::stoi(tmpInd);
            }

            std::stringstream matrixRowStringStream(tmpValues); 
            std::string matrixEntry;
            while (std::getline(matrixRowStringStream, matrixEntry, ' ')) {
                values.push_back(stod(matrixEntry));
            }

            // Assign values to robot model
            if (field == "numJoints") {
                numJoints = values[0]; 
                S.resize(numJoints);
                I.resize(numJoints);
                XTree.resize(numJoints);
                CoM.resize(numJoints);
                lam.resize(numJoints);
                transI.resize(numJoints);
                gearRatios.resize(numJoints);
                friction.resize(numJoints);
                damping.resize(numJoints);
            } else if (field == "twist") {
                S[ind] = Twist(values[0], values[1], values[2], values[3], values[4], values[5]);
            } else if (field == "gravity") {
                Eigen::Matrix<double, 3, 1> w, v;
                w << 0, 0, 0;
                v << values[0], values[1], values[2];
                gravity = Twist(w, v);
            } else if (field == "inertia") {
                I[ind].m = values[0];
                I[ind].I_bar << values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9]; 
                I[ind].m_c_hat << values[10], values[11], values[12], values[13], values[14], values[15], values[16], values[17], values[18];
            } else if (field == "Xtree") {
                XTree[ind].R << values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8];
                XTree[ind].p << values[9], values[10], values[11];
            } else if (field == "parent") {
                for (int j = 0; j < numJoints; j++) {
                    lam[j] = int(values[j]);
                }
            } else if (field == "CoM") {
                CoM[ind].p << values[0], values[1], values[2];
            } else if (field == "transI") {
                for (int j = 0; j < numJoints; j++) {
                    transI[j] = values[j];
                }
            } else if (field == "friction") {
                for (int j = 0; j < numJoints; j++) {
                    friction[j] = values[j];
                }
            } else if (field == "damping") {
                for (int j = 0; j < numJoints; j++) {
                    damping[j] = values[j];
                }
            } else if (field == "torque_limits") {

            } else if (field == "joint_limits") {

            } else if (field == "gear_ratios") {

            }
        }
    }

    std::vector<Twist> new_S; //Joint twists
    std::vector<RigidInertia> new_I; //Link spatial inertias
    std::vector<Transform> new_XTree; //Body tree transforms
    new_S.resize(numJoints);
    new_I.resize(numJoints);
    new_XTree.resize(numJoints);

    // Convert to Roy Featherstone style model
    for (int i = 0; i < numJoints; i++) {
        // Convert S from joint frame to world frame
        Transform Xwj = XTree[i];
        int pind = lam[i];
        
        while (pind > -1) {
            Xwj = Xwj.apply(XTree[pind]);
            pind = lam[pind];
        }
        new_S[i] = Xwj.invapply(S[i]);

        // Convert I from joint frame to body CoM frame
        new_I[i] = CoM[i].apply(I[i]);

        // Convert Xtree from joint-to-joint to CoM-to-CoM (backwards)
        Transform prevCoM;
        if (lam[i] != -1) {
            prevCoM = CoM[lam[i]];
        }

        Transform nextCoM = CoM[i];
        new_XTree[i] = prevCoM.apply(XTree[i].inverse().apply(nextCoM.inverse()));
    }

    S = new_S;
    I = new_I;
    XTree = new_XTree;
}

IntModel::IntModel(int numJoints_in) {
    numJoints = numJoints_in;

    S.resize(numJoints);
    I.resize(numJoints);
    XTree.resize(numJoints);
    CoM.resize(numJoints);
    lam.resize(numJoints);
    transI.resize(numJoints);
    gearRatios.resize(numJoints);
    friction.resize(numJoints);
    damping.resize(numJoints);
}

IntModel::IntModel(const Model* model, double eps) {
    numJoints = model->numJoints;

    S.resize(numJoints);
    I.resize(numJoints);
    XTree.resize(numJoints);
    CoM.resize(numJoints);
    lam.resize(numJoints);
    transI.resize(numJoints);
    gearRatios.resize(numJoints);
    friction.resize(numJoints);
    damping.resize(numJoints);

    // Populate interval vectors
    for (int i = 0; i < numJoints; i++) {
        // Copy S
        S[i].w = model->S[i].w.cast<Interval>();
        S[i].v = model->S[i].v.cast<Interval>();
        S[i].w_hat = model->S[i].w_hat.cast<Interval>();

        // Copy I
        I[i].m = Interval {model->I[i].m, model->I[i].m};
        I[i].I_bar = model->I[i].I_bar.cast<Interval>();
        I[i].m_c_hat = model->I[i].m_c_hat.cast<Interval>();

        // Copy XTree
        XTree[i].R = model->XTree[i].R.cast<Interval>();
        XTree[i].p = model->XTree[i].p.cast<Interval>();

        // Copy CoM
        CoM[i].p = model->CoM[i].p.cast<Interval>();
    }

    // Copy transI, gravity, gearRatios, friction, and damping
    transI = model->transI.cast<Interval>();
    gravity.w = model->gravity.w.cast<Interval>();
    gravity.v = model->gravity.v.cast<Interval>();
    gearRatios = model->gearRatios;
    friction = model->friction;
    damping = model->damping;

    // Copy parent array
    lam = model->lam;

    double lowP = 1 - eps;
    double highP = 1 + eps;
    
    // Add mass and inertia uncertainty
    for (int i = 0; i < numJoints; i++) {
        I[i].m.assign(I[i].m.lower()*lowP, I[i].m.lower()*highP); // Mass uncertainty

        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                double val = I[i].I_bar(j, k).lower();

                // Inertia uncertainty
                if (val >= 0) {
                    I[i].I_bar(j, k).assign(val*lowP, val*highP);
                }
                else { // If negative, need to flip to make sure lower < upper
                    I[i].I_bar(j, k).assign(val*highP, val*lowP);
                }
            }
        }
    }

    // Initialize RNEA storage variables
    v.resize(numJoints);
    va.resize(numJoints);
    a.resize(numJoints);
    f.resize(numJoints);
    Xbw.resize(numJoints);
    Sb.resize(numJoints);
    Xli_i.resize(numJoints);
}

Robot::Robot(std::string filename, double eps) {
    RobotModelPtr = new Model(filename);
    IntRobotModelPtr = new IntModel(RobotModelPtr, eps);
    numJoints = RobotModelPtr->numJoints;
}

Robot::~Robot() {
    delete RobotModelPtr;
    delete IntRobotModelPtr;
}

double Robot::integrateStateError(const Eigen::VectorXd &q, const Eigen::VectorXd &qd) {
    double stateError = 0;
    for (int i = 0; i < RobotModelPtr->numJoints; i++) {
        stateError += pow(qd(i) - q(i), 2);
    }
    return sqrt(stateError); 
};

Eigen::VectorXd Robot::finalizeTorque(Eigen::VectorXd& tau) {
    return tau;
}

// Virtual function for filling in unactuated joints
void Robot::fillUnactJoints(Eigen::VectorXd &q, Eigen::VectorXd &q_d, Eigen::VectorXd &q_dd) {}

#endif