#ifndef RNEA_CPP
#define RNEA_CPP

#include "rnea.hpp"

void passRNEA(Eigen::VectorXd& tau, Model* model, Eigen::VectorXd& q, Eigen::VectorXd& qd, 
              Eigen::VectorXd& qda, Eigen::VectorXd& qdd, bool applyFriction, bool applyGravity) {
    // Copy everything locally from model to improve readability
    int numJoints = model->numJoints;
    std::vector<Twist>& S = model->S; //Joint twists
    std::vector<RigidInertia>& I = model->I; //Link spatial inertias
    std::vector<Transform>& XTree = model->XTree; //Body tree transforms
    std::vector<int>& lam = model->lam; //Parent array that details joint connectivity
    Twist negGravity;

    if (applyGravity) {
        negGravity = -model->gravity;
    }

    // Declare velocity, acceleration and force arrays
    Twist v[numJoints];
    Twist va[numJoints];
    Twist a[numJoints];
    Wrench f[numJoints];

    // Storage variables
    Transform Xbw[numJoints]; //Transform from body to world frame
    Twist Sb[numJoints]; //Screw axis in body frame
    Transform Xli_i[numJoints]; //Transform from body i-1 to i

    // Forward pass
    for (int i = 0; i < numJoints; i++) {
        // Get parent index
        int li = lam[i];

        // Define ith body to world transform
        if (li != -1) {
            Xbw[i] = Xbw[lam[i]].apply(model->XTree[i]);
        }
        else {
            Xbw[i] = model->XTree[i];
        }

        // Screw axis in body frame
        Sb[i] = Xbw[i].invapply(S[i]);

        // Joint transform in body frame
        Xli_i[i] = Transform(Sb[i], -q[i]);

        // Transform from i - 1 to i
        Xli_i[i] = Xli_i[i].apply(model->XTree[i].inverse());

        if (li == -1) {
            // Velocity update from base
            v[i] = Sb[i]*qd[i];

            // Auxiliary velocity update from base
            va[i] = Sb[i]*qda[i];

            // Acceleration update from base
            a[i] = Xli_i[i].apply(negGravity) + Sb[i]*qdd[i] + 
                   v[i].cross(va[i]);
        }
        else {  
            // Velocity update
            v[i] = Xli_i[i].apply(v[li]) + Sb[i]*qd[i];

            // Auxiliary velocity update
            va[i] = Xli_i[i].apply(va[li]) + Sb[i]*qda[i];

            // Acceleration update
            a[i] = Xli_i[i].apply(a[li]) + Sb[i]*qdd[i] + 
                   v[i].cross(Sb[i]*qda[i]);
        }

        // Calculate v x Iv Jon's way
        Wrench vIv;
        vIv.tau = va[i].w.cross(I[i].I_bar * v[i].w);
        vIv.tau += I[i].I_bar * va[i].w.cross(v[i].w);
        vIv.f = I[i].m * va[i].w.cross(v[i].v);

        f[i] = I[i].apply(a[i]) + vIv;
    }

    // Backwards pass
    for (int i = numJoints - 1; i >= 0; i--) {
        tau[i] = Sb[i].dot(f[i]) + model->transI[i]*qdd[i];
        tau[i] += model->damping[i]*qd[i]; //Damping term
        if (applyFriction) tau[i] += model->friction[i]*( (qd[i] > 0) - (qd[i] < 0)); //This performs friction*signum(qd[i])
        if (lam[i] != -1) {
            f[lam[i]] = f[lam[i]] + Xli_i[i].invapply(f[i]);
        }
    }
}

void passRNEA_Int(VectorXint& tau, IntModel* model, Eigen::VectorXd& q, Eigen::VectorXd& qd, Eigen::VectorXd& qda, Eigen::VectorXd& qdd, bool applyFriction, bool applyGravity) {
    //Copy everything locally from model to improve readability
    int numJoints = model->numJoints;
    std::vector<IntTwist>& S = model->S; //Joint twists
    std::vector<IntRigidInertia>& I = model->I; //Link spatial inertias
    std::vector<IntTransform>& XTree = model->XTree; //Body tree transforms
    std::vector<int>& lam = model->lam; //Parent array that details joint connectivity
    IntTwist negGravity;

    //Storage variables for calculation
    std::vector<IntTwist>& v = model->v;
    std::vector<IntTwist>& va = model->va;
    std::vector<IntTwist>& a = model->a;
    std::vector<IntWrench>& f = model->f;
    std::vector<IntTransform>& Xbw = model->Xbw; //Transform from body to world frame
    std::vector<IntTwist>& Sb = model->Sb; //Screw axis in body frame
    std::vector<IntTransform>& Xli_i = model->Xli_i; //Transform from body i-1 to i


    if (applyGravity) {
        negGravity = -model->gravity;
    }

    //Forward pass
    for (int i = 0; i < numJoints; i++) {
        //Get parent index
        int li = lam[i];

        //Define ith body to world transform
        if (li != -1) {
            Xbw[i] = Xbw[lam[i]].apply(XTree[i]);
        }
        else {
            Xbw[i] = XTree[i];
        }

        //Screw axis in body frame
        Sb[i] = Xbw[i].invapply(S[i]);

        //Joint transform in body frame
        Xli_i[i] = IntTransform(Sb[i], -q[i]);

        //Transform from i - 1 to i
        Xli_i[i] = Xli_i[i].apply(XTree[i].inverse());

        IntTwist t;
        if (li == -1) {
            //Velocity update from base
            v[i] = Sb[i]*qd[i];

            //Auxiliary velocity update from base
            va[i] = Sb[i]*qda[i];

            //Acceleration update from base
            a[i] = Xli_i[i].apply(negGravity) + Sb[i]*qdd[i] + 
                   v[i].cross(va[i]);
        }
        else {  
            //Velocity update
            v[i] = Xli_i[i].apply(v[li]) + Sb[i]*qd[i];

            //Auxiliary velocity update
            IntTwist temp = Sb[i]*qda[i];
            va[i] = Xli_i[i].apply(va[li]) + temp;

            //Acceleration update
            a[i] = Xli_i[i].apply(a[li]) + Sb[i]*qdd[i] + 
                   v[i].cross(temp);
        }

        //Calculate v x Iv Jon's way
        IntWrench vIv;
        vIv.tau = va[i].w.cross(I[i].I_bar * v[i].w);
        vIv.tau += I[i].I_bar * va[i].w.cross(v[i].w);
        vIv.f = I[i].m * va[i].w.cross(v[i].v);

        f[i] = I[i].apply(a[i]) + vIv;
    }

    //Backwards pass
    for (int i = numJoints - 1; i >= 0; i--) {
        tau[i] = Sb[i].dot(f[i]);
        tau[i] += model->transI[i]*qdd[i]; // transmission inertia
        tau[i] += model->damping[i]*qd[i]; // damping
        if (applyFriction) tau[i] += model->friction[i]*( (qd[i] > 0) - (qd[i] < 0)); // static friction
        if (lam[i] != -1) {
            f[lam[i]] = f[lam[i]] + Xli_i[i].invapply(f[i]);
        }
    }
}

#endif
