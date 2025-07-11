#include "joint/Distance.h"
#include "rigidbody/RigidBody.h"

#include <Eigen/Dense>
#include <cstdio>

Distance::Distance() : Joint()
{

}

Distance::Distance(RigidBody* _body0, RigidBody* _body1, const Eigen::Vector3f& _r0, const Eigen::Vector3f& _r1, float _d) : Joint(_body0, _body1, _r0, Eigen::Quaternionf::Identity(), _r1, Eigen::Quaternionf::Identity(), kDistance), d(_d)
{
    dim = 1;
    J0.setZero(1, 6);
    J1.setZero(1, 6);
    J0Minv.setZero(1, 6);
    J1Minv.setZero(1, 6);
    phi.setZero(1);
    lambda.setZero(1);
}

void Distance::computeJacobian()
{
    // TODO Objective #2
    // Compute the Jacobians J0,J1, J0Minv, J1Min and constraint error phi
    // of the distance constraint.
    //

    Eigen::Vector3f x0 = body0->x + body0->q * r0;
    Eigen::Vector3f x1 = body1->x + body1->q * r1;

    Eigen::Vector3f dx = x1 - x0;

    if (dx.norm() < 1e-6f) {
        dx = Eigen::Vector3f(0,0,1); // Avoid division by zero
    }

    phi(0) = dx.norm() - d; //phi : distance between the two points minus the desired distance, is zero if the constraint is satisfied

    J0.block(0, 0, 1, 3) = -dx.transpose(); // linear part for body0
    J1.block(0, 0, 1, 3) = dx.transpose();  // Linear part for body1

    Eigen::Vector3f x0_cross = r0.cross(dx);
    Eigen::Vector3f x1_cross = r1.cross(dx);
    J0.block(0, 3, 1, 3) = x0_cross.transpose(); // angular part for body0
    J1.block(0, 3, 1, 3) = -x1_cross.transpose(); // angular part for body1

    J0Minv.block(0, 0, 1, 3) = J0.block(0, 0, 1, 3) * (1.0f / body0->mass);
    J0Minv.block(0, 3, 1, 3) = J0.block(0, 3, 1, 3) * body0->Iinv;

    J1Minv.block(0, 0, 1, 3) = J1.block(0, 0, 1, 3) * (1.0f / body1->mass);
    J1Minv.block(0, 3, 1, 3) = J1.block(0, 3, 1, 3) * body1->Iinv;
}