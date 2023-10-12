#ifndef RNEA_HPP
#define RNEA_HPP

#include "robot_models.hpp"

void passRNEA(Eigen::VectorXd& tau, Model* model, Eigen::VectorXd& q, Eigen::VectorXd& qd, Eigen::VectorXd& qda, Eigen::VectorXd& qdd, bool applyFriction, bool applyGravity = true);

void passRNEA_Int(VectorXint& tau, IntModel* model, Eigen::VectorXd& q, Eigen::VectorXd& qd, Eigen::VectorXd& qda, Eigen::VectorXd& qdd, bool applyFriction, bool applyGravity = true);

#endif