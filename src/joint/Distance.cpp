#include "joint/Distance.h"
#include "rigidbody/RigidBody.h"

#include <Eigen/Dense>

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

}

