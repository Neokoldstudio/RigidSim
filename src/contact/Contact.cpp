#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"

float Contact::mu = 0.8f;

static Eigen::Matrix3f hat(const Eigen::Vector3f& v)
{
    Eigen::Matrix3f vhat;
    vhat << 0, -v(2), v(1),
        v(2), 0, -v(0),
        -v(1), v(0), 0;
    return vhat;
}

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
    Eigen::Vector3f t = n.cross(Eigen::Vector3f(0.0f, 0.0f, 1.0f));
    if (t.norm() < 1e-6f)
    {
        t = n.cross(Eigen::Vector3f(1.0f, 0.0f, 0.0f));
    }
    t.normalize();
    // TODO Compute second tangent direction b.
    //
    Eigen::Vector3f b = n.cross(t);
    b.normalize();
    // TODO Compute the Jacobians blocks J0 and J1
    //
    Eigen::Matrix3f ntb;
    ntb.col(0) = n;
    ntb.col(1) = t;
    ntb.col(2) = b;
    
    Eigen::Vector3f rr0 = body0->q * (p - body0->x);//r0 vector from body0's center of mass to the contact point transformed in world space
    Eigen::Vector3f rr1 = body1->q * (p - body1->x);//same for body1

    //Here we assemble the Jacobian blocks J0 and J1, both 3x6 matrices. The complete Jacobian matrix J is a 3x12 matrix as seen in the course slides.
    J0.block(0, 0, 3, 3) = -ntb; // Linear part for body0
    J0.block(0, 3, 3, 3) = hat(rr0) * ntb; // Angular part for body0
    J1.block(0, 0, 3, 3) = ntb; //same
    J1.block(0, 3, 3, 3) = -hat(rr1) * ntb; //same

    // TODO Finally, compute the blocks J M^-1 for each body.
    //
    J0Minv.block(0, 0, 3, 3) = (1.0f / body0->mass) * J0.block(0, 0, 3, 3);
    J0Minv.block(0, 3, 3, 3) = J0.block(0, 3, 3, 3) * body0->Iinv;
    J1Minv.block(0, 0, 3, 3) = (1.0f / body1->mass) * J1.block(0, 0, 3, 3);
    J1Minv.block(0, 3, 3, 3) = J1.block(0, 3, 3, 3) * body1->Iinv;
}
