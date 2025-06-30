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

    phi(0) = (x1 - x0).norm() - d; //phi : distance between the two points minus the desired distance, is zero if the constraint is satisfied

    J0.block(0, 0, 1, 3) = (x0 - x1).transpose(); // Distance Jacobian for body0 is the transposed vector from body0 to body1
    J1.block(0, 0, 4, 6) = (x1 - x0).transpose(); // Same, but the other way around
    
    J0Minv.block(0, 0, 1, 3) = J0.block(0, 0, 1, 3) * (1 / body0->mass);
    J0Minv.block(0, 0, 4, 6) = J0.block(0, 0, 4, 6) * body0->Iinv;

    J1Minv.block(0, 0, 1, 3) = J1.block(0, 0, 1, 3) * (1 / body1->mass);
    J1Minv.block(0, 0, 4, 6) = J1.block(0, 0, 4, 6) * body1->Iinv;
}