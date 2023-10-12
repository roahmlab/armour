#ifndef SPATIAL
#define SPATIAL

#include "headers.hpp"

class Wrench {
public:
    Eigen::Matrix<double, 3, 1> tau;
    Eigen::Matrix<double, 3, 1> f;

    Wrench(Eigen::Matrix<double, 3, 1> tau_in, Eigen::Matrix<double, 3, 1> f_in) {
        tau = tau_in;
        f = f_in;
    }

    Wrench(double t1, double t2, double t3, double f1, double f2, double f3) {
        tau << t1, t2, t3;
        f << f1, f2, f3;
    }
    
    Wrench() {
        tau = Eigen::Matrix<double, 3, 1>::Zero();
        f = Eigen::Matrix<double, 3, 1>::Zero();
    }

    friend std::ostream& operator<<(std::ostream& os, const Wrench& wrench) {
        os << wrench.tau << std::endl << wrench.f;
        return os;
    }

    Wrench operator+(const Wrench& w2) {
        Eigen::Matrix<double, 3, 1> newTau = tau + w2.tau;
        Eigen::Matrix<double, 3, 1> newF = f + w2.f;
        return Wrench(newTau, newF);
    }
};

class Twist {
public:
    Eigen::Matrix<double, 3, 1> w;
    Eigen::Matrix<double, 3, 1> v;
    Eigen::Matrix<double, 3, 3> w_hat;

    Twist(Eigen::Matrix<double, 3, 1> w_in, Eigen::Matrix<double, 3, 1> v_in) {
        w = w_in;
        v = v_in;
        w_hat << 0, -w(2), w(1),
                  w(2), 0, -w(0),
                  -w(1), w(0), 0;
    }

    Twist(double w1, double w2, double w3, double q1, double q2, double q3) {
        Eigen::Matrix<double, 3, 1> w_in{w1, w2, w3};
        Eigen::Matrix<double, 3, 1> q_in{q1, q2, q3};
        // define(w_in, q_in);
        w = w_in;
        v = q_in;
        w_hat << 0, -w(2), w(1),
                  w(2), 0, -w(0),
                  -w(1), w(0), 0;
    }

    Twist() {
        w = Eigen::Matrix<double, 3, 1>::Zero();
        v = Eigen::Matrix<double, 3, 1>::Zero();
        w_hat = Eigen::Matrix<double, 3, 3>::Zero();
    }

    void define(Eigen::Matrix<double, 3, 1> w_in, Eigen::Matrix<double, 3, 1> q_in) {
        w = w_in;
        w.normalize();
        v = -w.cross(q_in);
        w_hat << 0, -w(2), w(1),
                 w(2), 0, -w(0),
                 -w(1), w(0), 0;
    }

    double dot(Wrench& f) {
        double result = w.adjoint() * f.tau;
        return result + v.adjoint() * f.f;
    }

    Twist cross(Twist z2) {
        Eigen::Matrix<double, 3, 1> newW = w_hat * z2.w;
        Eigen::Matrix<double, 3, 1> newV = w_hat * z2.v + v.cross(z2.w);
        return Twist(newW, newV);
    }

    Wrench cross(Wrench w) {
        Eigen::Matrix<double, 3, 1> newTau = w_hat * w.tau + v.cross(w.f);
        Eigen::Matrix<double, 3, 1> newf = w_hat * w.f;
        return Wrench(newTau, newf);
    }

    friend std::ostream& operator<<(std::ostream& os, const Twist& zeta) {
        os << zeta.w << std::endl << zeta.v;
        return os;
    }

    Twist operator*(const double s) {
        Eigen::Matrix<double, 3, 1> newW = w*s;
        Eigen::Matrix<double, 3, 1> newV = v*s;
        return Twist(newW, newV);
    }

    Twist operator+(const Twist& z2) {
        Eigen::Matrix<double, 3, 1> newW = w + z2.w;
        Eigen::Matrix<double, 3, 1> newV = v + z2.v;
        return Twist(newW, newV);
    }

    Twist operator-() {
        return Twist(-w, -v);
    }
};

class RigidInertia {
public:
    double m; //Mass
    Eigen::Matrix<double, 3, 3> I_bar; //Inertia (top-left block)
    Eigen::Matrix<double, 3, 3> m_c_hat; //m*cx (bottom-left, top-right block)

    RigidInertia(double m_in, Eigen::Matrix<double, 3, 1> c, Eigen::Matrix<double, 3, 3> Ic) {
        m = m_in;
        Eigen::Matrix<double, 3, 3> c_hat;
        c_hat << 0, -c(2), c(1),
                 c(2), 0, -c(0),
                 -c(1), c(0), 0;
        m_c_hat = m*c_hat;
        I_bar = Ic - m_c_hat*c_hat;
    }

    RigidInertia(Eigen::Matrix<double, 3, 3> I_bar_in, Eigen::Matrix<double, 3, 3> m_c_hat_in) {
        I_bar = I_bar_in;
        m_c_hat = m_c_hat_in;
    }

    RigidInertia() {
        I_bar = Eigen::Matrix<double, 3, 3>::Zero();
        m_c_hat = Eigen::Matrix<double, 3, 3>::Zero();
    }

    Wrench apply(Twist zeta) {
        Eigen::Matrix<double, 3, 1> tau = I_bar*zeta.w + m_c_hat*zeta.v;
        Eigen::Matrix<double, 3, 1> f = m*zeta.v - m_c_hat*zeta.w;
        return Wrench(tau, f);
    }
};

class Transform {
public:
    Eigen::Matrix<double, 3, 3> R; //Rotation matrix
    Eigen::Matrix<double, 3, 1> p; //Translation vector
    Eigen::Matrix<double, 3, 3> I; //3x3 identity

    Transform(Twist zeta, double theta) {
        //Initialize identity
        I = Eigen::Matrix<double, 3, 3>::Identity();

        /* Define rotation matrix using w and theta */
        //Calculate R using Rodrigues' formula
        R = Eigen::Matrix<double, 3, 3>::Identity() + 
                                        zeta.w_hat*std::sin(theta) + 
                                        (1 - std::cos(theta)) * zeta.w_hat * zeta.w_hat;

        /* Calculate p using (I - R) * (w x v). */
        //Leaving out the ww'v*theta term because v and w are orthogonal for revolute joints
        p = (I - R) * zeta.w_hat * zeta.v;
        p = -R.transpose()*p; //Featherstone convention
    }

    Transform(Eigen::Matrix<double, 3, 3> R_in, Eigen::Matrix<double, 3, 1> p_in) {
        //Initialize identity
        I = Eigen::Matrix<double, 3, 3>::Identity();
        R = R_in;
        p = p_in;
    }

    Transform() {
        //Initialize identity
        I = Eigen::Matrix<double, 3, 3>::Identity();
        R = I;
        p = Eigen::Matrix<double, 3, 1>::Zero();
    }

    //Apply transform using adjoints, v1 = Xv2, defined in
    //Effi.pdf of Featherstone's lecture notes
    Twist apply(Twist zeta) {
        Eigen::Matrix<double, 3, 1> newW = R*zeta.w;
        Eigen::Matrix<double, 3, 1> newV = R*(zeta.v - p.cross(zeta.w));
        return Twist(newW, newV);
    }

    //Apply the inverse of the transform, defined in
    //Effi.pdf of Featherstone's lecture notes
    Twist invapply(Twist& zeta) {
        Eigen::Matrix<double, 3, 1> newW = R.transpose()*zeta.w;
        Eigen::Matrix<double, 3, 1> newV = R.transpose()*zeta.v + p.cross(newW);
        return Twist(newW, newV);
    }

    //Apply wrench transform, f1 = X*f2, where X* = X^(-T) (inverse transpose)
    //
    Wrench apply(Wrench& w) {
        Eigen::Matrix<double, 3, 1> newTau = R*(w.tau - p.cross(w.f));
        Eigen::Matrix<double, 3, 1> newF = R*w.f;
        return Wrench(newTau, newF);
    }

    Wrench invapply(Wrench& w) {
        Eigen::Matrix<double, 3, 1> newTau = R.transpose()*w.tau + p.cross(R.transpose()*w.f);
        Eigen::Matrix<double, 3, 1> newF = R.transpose()*w.f;
        return Wrench(newTau, newF);
    }

    Eigen::Matrix<double, 3, 1> apply(Eigen::Matrix<double, 3, 1> v) {
        return R*v + p;
    }

    RigidInertia apply(RigidInertia I) {
        Eigen::Matrix<double, 3, 3> p_hat;
        p_hat << 0, -p(2), p(1),
                 p(2), 0, -p(0),
                 -p(1), p(0), 0;

        RigidInertia newI;
        Eigen::Matrix<double, 3, 3> mRp_hat = I.m*R*p_hat;

        newI.m = I.m;
        // R*h*R' - m*R*C_hat*R';
        newI.m_c_hat = R*I.m_c_hat*R.transpose() - I.m*R*p_hat*R.transpose();
        newI.I_bar = (R*(I.I_bar + 2*I.m_c_hat*p_hat) - mRp_hat*p_hat)*R.transpose();
        return newI;
    }

    Transform apply(Transform x2) {
        //This definition of X1*X2 works for twists but not for
        //homogeneous position vectors. Not sure why.
        Eigen::Matrix<double, 3, 3> newR = R*x2.R;
        Eigen::Matrix<double, 3, 1> newP = x2.p + x2.R.transpose()*p;

        return Transform(newR, newP);
    }

    Transform inverse() {
        Transform newT;
        newT.R = R.transpose();
        newT.p = -R*p;  
        return newT;
    }
};

#endif
