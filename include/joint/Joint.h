#pragma once

#include <Eigen/Dense>

#include "util/Types.h"

class RigidBody;

enum eJointType { kNoJoint = -1, kSpherical = 0, kDistance };

// Joint class.
//
class Joint
{
public:

    // Constructor with all parameters.
    Joint(RigidBody* _body0, RigidBody* _body1, const Eigen::Vector3f& _r0, const Eigen::Quaternionf& _q0, const Eigen::Vector3f& _r1, const Eigen::Quaternionf& _q1, eJointType _type);

    virtual ~Joint() { }

    RigidBody* body0;           // The first body
    RigidBody* body1;           // The second body
    JBlock J0;                  // The constraint Jacobian of body0
    JBlock J1;                  // The constraint Jacobian of body1
    JBlock J0Minv;              // The constraint Jacobian of body0 * inverse Mass
    JBlock J1Minv;              // The constraint Jacobian of body1 * inverse Mass
    Eigen::VectorXf phi;        // Contraint error
    Eigen::VectorXf lambda;     // Constraint impulse


    // The constraint attachment.
    // The global attachment point on body0 can be computed by:
    // 
    //    r0_global = (body0->q * r0) + body0->x
    // 
    // and the orientation at the attachment is simply:
    //    q0_global = body0->q * q0
    //
    // The steps are similar for the second attachment, but using body1.
    //
    Eigen::Vector3f r0;         // Relative attachment point of joint in body0 coordinate frame.
    Eigen::Vector3f r1;         // Relative attachment point of joint in body1 coordinate frame.
    Eigen::Quaternionf q0;      // Relative attachment orientation in body0 coordinate frame.
    Eigen::Quaternionf q1;      // Relative attachment orientation in body1 coordinate frame.
    
    unsigned int dim;           // Number of constraint equations, aka # of degrees of freedom (dof) removed 
            
    eJointType type;            // The joint type.

    virtual void computeJacobian();

protected:

    // Default constructor (hidden).
    Joint();

};
