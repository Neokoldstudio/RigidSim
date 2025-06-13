#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"

float Contact::mu = 0.8f;

Contact::Contact() : p(), n(), t(), b()
{

}

Contact::Contact(RigidBody* _body0, RigidBody* _body1, const Eigen::Vector3f& _p, const Eigen::Vector3f& _n, float _pene) :
    body0(_body0), body1(_body1),
    p(_p), n(_n), t(), b(), pene(_pene)
{
    J0.setZero();
    J1.setZero();
    J0Minv.setZero();
    J1Minv.setZero();
    lambda.setZero();
    phi.setZero();
    phi(0) = _pene;

    // Track list of contacts for each body.
    // This makes it easy to find coupled constraints and contacts during the solve.
    //
    // Note: Unlike joints, which have the list maintained by RigidBodySystem::addJoint()
    //  here we update the list in the constructor since contacts are determined during simulation.
    //
    body0->contacts.push_back(this);
    body1->contacts.push_back(this);
}

Contact::~Contact()
{

}

void Contact::computeJacobian()
{
    // TODO Objective #5
    // Compute the Jacobian for a contact constraint.
    // 
   
    // TODO Compute the contact frame, 
    // which consists of an orthonormal bases formed the vector n, t, and b
    //
    // Compute first tangent direction t
    //

    // TODO Compute second tangent direction b.
    //

    // TODO Compute the Jacobians blocks J0 and J1
    //

    // TODO Finally, compute the blocks J M^-1 for each body.
    //
}
