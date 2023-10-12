#ifndef SPATIAL_INTERVAL
#define SPATIAL_INTERVAL

#include "headers.hpp"

class IntWrench {
public:
    Eigen::Matrix<Interval, 3, 1> tau;
    Eigen::Matrix<Interval, 3, 1> f;

    IntWrench(Eigen::Matrix<Interval, 3, 1> tau_in, Eigen::Matrix<Interval, 3, 1> f_in) {
        tau = tau_in;
        f = f_in;
    }

    IntWrench(Interval t1, Interval t2, Interval t3, Interval f1, Interval f2, Interval f3) {
        tau << t1, t2, t3;
        f << f1, f2, f3;
    }
    
    IntWrench() {
        tau = Eigen::Matrix<Interval, 3, 1>::Zero();
        f = Eigen::Matrix<Interval, 3, 1>::Zero();
    }

    IntWrench operator+(const IntWrench& w2) {
        Eigen::Matrix<Interval, 3, 1> newTau = tau + w2.tau;
        Eigen::Matrix<Interval, 3, 1> newF = f + w2.f;
        return IntWrench(newTau, newF);
    }
};

class IntTwist {
public:
    Eigen::Matrix<Interval, 3, 1> w;
    Eigen::Matrix<Interval, 3, 1> v;
    Eigen::Matrix<Interval, 3, 3> w_hat;

    IntTwist(Eigen::Matrix<Interval, 3, 1> &w_in, Eigen::Matrix<Interval, 3, 1> &v_in) {
        w = w_in;
        v = v_in;
        w_hat << 0, -w(2), w(1),
                  w(2), 0, -w(0),
                  -w(1), w(0), 0;
    }

    IntTwist(Interval &w1, Interval &w2, Interval &w3, Interval &q1, Interval &q2, Interval &q3) {
        Eigen::Matrix<Interval, 3, 1> w{w1, w2, w3};
        Eigen::Matrix<Interval, 3, 1> q{q1, q2, q3};
        define(w, q);
    }

    IntTwist() {
        w = Eigen::Matrix<Interval, 3, 1>::Zero();
        v = Eigen::Matrix<Interval, 3, 1>::Zero();
        w_hat = Eigen::Matrix<Interval, 3, 3>::Zero();
    }

    void define(Eigen::Matrix<Interval, 3, 1> &w_in, Eigen::Matrix<Interval, 3, 1> &q_in) {
        w = w_in;
        w.normalize();
        v = -w.cross(q_in);
        w_hat << 0, -w(2), w(1),
                 w(2), 0, -w(0),
                 -w(1), w(0), 0;
    }

    Interval dot(IntWrench& f) {
        Interval result = w.adjoint() * f.tau;
        Interval result2 = v.adjoint() * f.f;
        return result + result2;
    }

    IntTwist cross(IntTwist &z2) {
        Eigen::Matrix<Interval, 3, 1> newW = w_hat * z2.w;
        Eigen::Matrix<Interval, 3, 1> newV = w_hat * z2.v + v.cross(z2.w);
        return IntTwist(newW, newV);
    }

    IntWrench cross(IntWrench &w) {
        Eigen::Matrix<Interval, 3, 1> newTau = w_hat * w.tau + v.cross(w.f);
        Eigen::Matrix<Interval, 3, 1> newf = w_hat * w.f;
        return IntWrench(newTau, newf);
    }

    IntTwist operator*(const Interval &s) {
        Eigen::Matrix<Interval, 3, 1> newW = w*s;
        Eigen::Matrix<Interval, 3, 1> newV = v*s;
        return IntTwist(newW, newV);
    }

    IntTwist operator+(const IntTwist& z2) {
        Eigen::Matrix<Interval, 3, 1> newW = w + z2.w;
        Eigen::Matrix<Interval, 3, 1> newV = v + z2.v;
        return IntTwist(newW, newV);
    }

    IntTwist operator-() {
        Eigen::Matrix<Interval, 3, 1> neg_w = -w;
        Eigen::Matrix<Interval, 3, 1> neg_v = -v;
        return IntTwist(neg_w, neg_v);
    }
};

class IntRigidInertia {
public:
    Interval m; //Mass
    Eigen::Matrix<Interval, 3, 3> I_bar; //Inertia (top-left block)
    Eigen::Matrix<Interval, 3, 3> m_c_hat; //m*cx (bottom-left, top-right block)

    IntRigidInertia(Interval &m_in, Eigen::Matrix<Interval, 3, 1> &c, Eigen::Matrix<Interval, 3, 3> &Ic) {
        m = m_in;
        Eigen::Matrix<Interval, 3, 3> c_hat;
        c_hat << 0, -c(2), c(1),
                 c(2), 0, -c(0),
                 -c(1), c(0), 0;
        m_c_hat = m*c_hat;
        I_bar = Ic - m_c_hat*c_hat;
    }

    IntRigidInertia(Eigen::Matrix<Interval, 3, 3> &I_bar_in, Eigen::Matrix<Interval, 3, 3> &m_c_hat_in) {
        I_bar = I_bar_in;
        m_c_hat = m_c_hat_in;
    }

    IntRigidInertia() {
        I_bar = Eigen::Matrix<Interval, 3, 3>::Zero();
        m_c_hat = Eigen::Matrix<Interval, 3, 3>::Zero();
    }

    IntWrench apply(const IntTwist &zeta) {
        Eigen::Matrix<Interval, 3, 1> tau = I_bar*zeta.w + m_c_hat*zeta.v;
        Eigen::Matrix<Interval, 3, 1> f = m*zeta.v - m_c_hat*zeta.w;
        return IntWrench(tau, f);
    }
};

const Eigen::Matrix<Interval, 3, 3> interval_identity = Eigen::Matrix<Interval, 3, 3>::Identity();
class IntTransform {
public:
    Eigen::Matrix<Interval, 3, 3> R; //Rotation matrix
    Eigen::Matrix<Interval, 3, 1> p; //Translation vector
    static Eigen::Matrix<Interval, 3, 3> I;

    IntTransform(IntTwist &zeta, double theta) {
        /* Define rotation matrix using w and theta */
        //Calculate R using Rodrigues' formula
        R = interval_identity + zeta.w_hat*std::sin(theta) + 
              (1 - std::cos(theta)) * zeta.w_hat * zeta.w_hat;

        /* Calculate p using (I - R) * (w x v). */
        //Leaving out the ww'v*theta term because v and w are orthogonal for revolute joints
        p = (interval_identity - R) * zeta.w_hat * zeta.v;
        p = -R.transpose()*p; //Featherstone convention
    }

    IntTransform(Eigen::Matrix<Interval, 3, 3> &R_in, Eigen::Matrix<Interval, 3, 1> &p_in) {
        R = R_in;
        p = p_in;
    }

    IntTransform() {
        R = interval_identity;
        p = Eigen::Matrix<Interval, 3, 1>::Zero();
    }

    //Apply transform using adjoints, v1 = Xv2, defined in
    //Effi.pdf of Featherstone's lecture notes
    IntTwist apply(const IntTwist &zeta) {
        Eigen::Matrix<Interval, 3, 1> newW = R*zeta.w;
        Eigen::Matrix<Interval, 3, 1> newV = R*(zeta.v - p.cross(zeta.w));
        return IntTwist(newW, newV);
    }

    //Apply the inverse of the transform, defined in
    //Effi.pdf of Featherstone's lecture notes
    IntTwist invapply(const IntTwist& zeta) {
        Eigen::Matrix<Interval, 3, 1> newW = R.transpose()*zeta.w;
        Eigen::Matrix<Interval, 3, 1> newV = R.transpose()*zeta.v + p.cross(newW);
        return IntTwist(newW, newV);
    }

    //Apply IntWrench transform, f1 = X*f2, where X* = X^(-T) (inverse transpose)
    //
    IntWrench apply(const IntWrench& w) {
        Eigen::Matrix<Interval, 3, 1> newTau = R*(w.tau - p.cross(w.f));
        Eigen::Matrix<Interval, 3, 1> newF = R*w.f;
        return IntWrench(newTau, newF);
    }

    IntWrench invapply(const IntWrench& w) {
        Eigen::Matrix<Interval, 3, 1> newTau = R.transpose()*w.tau + p.cross(R.transpose()*w.f);
        Eigen::Matrix<Interval, 3, 1> newF = R.transpose()*w.f;
        return IntWrench(newTau, newF);
    }

    Eigen::Matrix<Interval, 3, 1> apply(const Eigen::Matrix<Interval, 3, 1> &v) {
        return R*v + p;
    }

    IntRigidInertia apply(const IntRigidInertia &inertia) {
        Eigen::Matrix<Interval, 3, 3> p_hat;
        p_hat << 0, -p(2), p(1),
                 p(2), 0, -p(0),
                 -p(1), p(0), 0;

        IntRigidInertia newI;
        Eigen::Matrix<Interval, 3, 3> mRp_hat = inertia.m*R*p_hat;

        newI.m = inertia.m;
        // R*h*R' - m*R*C_hat*R';
        newI.m_c_hat = R*inertia.m_c_hat*R.transpose() - inertia.m*R*p_hat*R.transpose();
        newI.I_bar = (R*(inertia.I_bar + 2*inertia.m_c_hat*p_hat) - mRp_hat*p_hat)*R.transpose();
        return newI;
    }

    IntTransform apply(const IntTransform &x2) {
        //This definition of X1*X2 works for twists but not for
        //homogeneous position vectors. Not sure why.
        Eigen::Matrix<Interval, 3, 3> newR = R*x2.R;
        Eigen::Matrix<Interval, 3, 1> newP = x2.p + x2.R.transpose()*p;

        return IntTransform(newR, newP);
    }

    IntTransform inverse() {
        IntTransform newT;
        newT.R = R.transpose();
        newT.p = -R*p;  
        return newT;
    }
};

#endif
