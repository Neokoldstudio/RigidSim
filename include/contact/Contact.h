#pragma once

#include "joint/Joint.h"

#include <Eigen/Dense>

class RigidBody;

// Contact constraint with box friction.
// 
// This class stores the contact normal @a n and contact point @a p,
// as well as the Jacobians @a J0 and @a J1 for each body,
// which are computed by computeJacobian().
// 
class Contact
{
public: 
    static float mu;            // Coefficient of friction (global)

public:

    // Constructor with all parameters.
    Contact(RigidBody* _body0, RigidBody* _body1, const Eigen::Vector3f& _p, const Eigen::Vector3f& _n, float _pene);

    virtual ~Contact();

    RigidBody* body0;
    RigidBody* body1;
    JCBlock J0, J1;
    JCBlock J0Minv, J1Minv;
    Eigen::Vector3f lambda, phi;

    Eigen::Vector3f p;          // The contact point.
    Eigen::Vector3f n;          // The contact normal.
    Eigen::Vector3f t, b;       // Tangent directions.
    float pene;                 // Penetration

    // Computes the constraint Jacobian used for contact with friction. 
    // 
    // This function also builds the contact frame using the contact normal @a n
    // and the provided direction @a dir, which is aligned with the first tangent direction.
    // The resulting frame is stored in the basis vectors @a n, @a t1, and @a t2.
    //
    void computeJacobian();


protected:

    // Default constructor (hidden)
    explicit Contact();
};
