#include "joint/Joint.h"
#include "rigidbody/RigidBody.h"

#include <memory>

namespace
{
    static const Eigen::Matrix3f sEye = Eigen::Matrix3f::Identity();

    static inline Eigen::Matrix3f hat(const Eigen::Vector3f& v)
    {
        Eigen::Matrix3f vhat;
        vhat << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0), 0;
        return vhat;
    }

    static inline void computeHingeJacobian(Joint* j)
    {
        const Eigen::Vector3f rr0 = j->body0->q * j->r0;
        const Eigen::Vector3f rr1 = j->body1->q * j->r1;
        const Eigen::Vector3f nn = j->body0->q * (j->q0 * Eigen::Vector3f(1, 0, 0));
        const Eigen::Vector3f uu = j->body1->q * (j->q1 * Eigen::Vector3f(0, 1, 0));
        const Eigen::Vector3f vv = j->body1->q * (j->q1 * Eigen::Vector3f(0, 0, 1));
        const Eigen::Vector3f ncrossuu = nn.cross(uu);
        const Eigen::Vector3f ncrossvv = nn.cross(vv);

        // compute constraint error
        j->phi.head(3) = (j->body0->x + rr0 - j->body1->x - rr1);
        j->phi(3) = nn.dot(uu);
        j->phi(4) = nn.dot(vv);

        // Compute Jacobian
        j->J0.block(0, 0, 3, 3) = sEye;
        j->J0.block(3, 3, 1, 3) = ncrossuu.transpose();
        j->J0.block(4, 3, 1, 3) = ncrossvv.transpose();
        j->J0.block(0, 3, 3, 3) = hat(-rr0);
        j->J1.block(0, 0, 3, 3) = -sEye;
        j->J1.block(0, 3, 3, 3) = hat(rr1);
        j->J1.block(3, 3, 1, 3) = -ncrossuu.transpose();
        j->J1.block(4, 3, 1, 3) = -ncrossvv.transpose();

        j->J0Minv.block(0, 0, 5, 3) = (1.0f / j->body0->mass) * j->J0.block(0, 0, 5, 3);
        j->J0Minv.block(0, 3, 5, 3) = j->J0.block(0, 3, 5, 3) * j->body0->Iinv;
        j->J1Minv.block(0, 0, 5, 3) = (1.0f / j->body1->mass) * j->J1.block(0, 0, 5, 3);
        j->J1Minv.block(0, 3, 5, 3) = j->J1.block(0, 3, 5, 3) * j->body1->Iinv;
    }

    static inline void computeSphericalJacobian(Joint* j)
    {
        const Eigen::Vector3f rr0 = j->body0->q * j->r0;
        const Eigen::Vector3f rr1 = j->body1->q * j->r1;

        // compute constraint error
        j->phi.head(3) = (j->body0->x + rr0 - j->body1->x - rr1);

        // Compute Jacobian
        j->J0.block(0, 0, 3, 3) = sEye;
        j->J1.block(0, 0, 3, 3) = -sEye;
        j->J0.block(0, 3, 3, 3) = hat(-rr0);
        j->J1.block(0, 3, 3, 3) = hat(rr1);
        j->J0Minv.block(0, 0, 3, 3) = (1.0f / j->body0->mass) * j->J0.block(0, 0, 3, 3);
        j->J0Minv.block(0, 3, 3, 3) = j->J0.block(0, 3, 3, 3) * j->body0->Iinv;
        j->J1Minv.block(0, 0, 3, 3) = (1.0f / j->body1->mass) * j->J1.block(0, 0, 3, 3);
        j->J1Minv.block(0, 3, 3, 3) = j->J1.block(0, 3, 3, 3) * j->body1->Iinv;
    }

    static inline void noop(Joint* j)
    {

    }

}


Joint::Joint() :
    body0(nullptr), body1(nullptr), dim(0), type(kNoJoint)
{

}

Joint::Joint(RigidBody* _body0, RigidBody* _body1, 
    const Eigen::Vector3f& _r0, const Eigen::Quaternionf& _q0, const Eigen::Vector3f& _r1, const Eigen::Quaternionf& _q1, eJointType _type) :
    body0(_body0), body1(_body1), dim(0), type(_type), r0(_r0), r1(_r1), q0(_q0), q1(_q1)
{
    J0.setZero();
    J1.setZero();
    J0Minv.setZero();
    J1Minv.setZero();
    lambda.setZero();
    phi.setZero();
}

void Joint::computeJacobian()
{

}
