#include "solvers/SolverPGS.h"

#include "contact/Contact.h"
#include "joint/Joint.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

#include <Eigen/Dense>


namespace
{
    // TODO 
    // Solve the Boxed LCP problem for a single contact and isotropic Coulomb friction.
    // The solution vector, @a x, contains the impulse the non-interpenetration constraint in x(0), and
    // the friction constraints in x(1) and x(2)
    // 
    // The solution of the non-interpenetration impulses must be positive.
    // The solution of friction must projected to the lower and upper bounds imposed by the box model.
    // 
    // Inputs: 
    //    A - diagonal block of the global JMinvJT matrix corresponding to this contact
    //    x - contains three impulse variables (non-interpenetration + two friction)
    //    b - updated rhs vector minus contributions of coupled joints and contacts
    //    mu - the friction coefficient
    //
    static inline void solveContact(const Eigen::Matrix3f& A, const Eigen::Vector3f& b, Eigen::Vector3f& x, const float mu)
    {

    }

    // TODO 
    // Solve for the constraint impluses of a bilateral constraint.
    // The solution vector, @a x, contains the impulse the non-interpenetration constraint in x(0), and
    // the friction constraints in x(1) and x(2)
    // 
    // The solution of the non-interpenetration impulses must be positive.
    // The solution of friction must projected to the lower and upper bounds imposed by the box model.
    // 
    // Inputs: 
    //    A - diagonal block of the global JMinvJT matrix corresponding to this joint
    //    x - contains the constraint impulses (lambda)
    //    b - updated rhs vector minus contributions of coupled joints and contacts
    //
    template <typename VectorBlockType>
    static inline void solveJoint(const Eigen::MatrixXf& A, const Eigen::VectorXf& b, VectorBlockType& x)
    {

    }
}

SolverPGS::SolverPGS(RigidBodySystem* _rigidBodySystem) : Solver(_rigidBodySystem)
{

}

void SolverPGS::solve(float h)
{
    std::vector<Contact>& contacts = m_rigidBodySystem->getContacts();
    std::vector<Joint*>& joints = m_rigidBodySystem->getJoints();
    const int numContacts = contacts.size();
    const int numJoints = joints.size();

    // TODO Objetive #6
    // Implement the matrix-free PGS algorithm to solve for constraint and contact impulses.
    // 
    // See TODO notes below, and the helper functions in the anonymous namespace above.
    //

    if (numJoints > 0)
    {
        // TODO Build diagonal matrices of bilateral joints
        //  such that each Ajoint[i] = J_i * Minv * tanspose(J_i)
        //

        // Build diagonal matrices
        Ajoint.resize(numJoints);
        for (int i = 0; i < numJoints; ++i)
        {
            Joint* j = joints[i];
            const int dim = j->dim;
            const float eps = 1e-6f;

            // Compute the diagonal term : Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            //
            Eigen::MatrixXf A = eps * Eigen::MatrixXf::Identity(dim, dim);

            if (!j->body0->fixed)
            {
                A += j->J0Minv.block(0, 0, dim, 6) * j->J0.block(0, 0, dim, 6).transpose();
            }
            if (!j->body1->fixed)
            {
                A += j->J1Minv.block(0, 0, dim, 6) * j->J1.block(0, 0, dim, 6).transpose();
            }

            Ajoint[i] = A;
        }
    }

    if (numContacts > 0)
    {
        Acontact.resize(numContacts);
        for (int i = 0; i < numContacts; ++i)
        {
            // TODO Build array of 3x3 diagonal matrices, one for each contact. 
            // Store the result in Acontact[i], such that each Acontact[i] = J_i * Minv * tanspose(J_i)
            // 

        }
    }

    bcontact.resize(numContacts);
    for (int i = 0; i < numContacts; ++i)
    {
        Contact* c = &contacts[i];
        c->lambda.setZero();

        // TODO Compute the right-hand side vector for contacts, e.g.
        //      b = -gamma*phi/h - J*vel - h*JMinvJT*force
        //
    }

    bjoint.resize(numJoints);
    for (int i = 0; i < numJoints; ++i)
    {
        Joint* j = joints[i];

        // TODO Compute the right-hand side vector for joints, e.g.
        //      b = -gamma*phi/h - J*vel - h*JMinvJT*force
        //

    }

    // TODO 
    // PGS main loop.
    // There is no convergence test here. Simply stop after @a maxIter iterations.
    //
    for(int iter = 0; iter < m_maxIter; ++iter)
    {
        // TODO
        // For each joint, compute an updated value of joints[i]->lambda.
        // 
        // IMPORTANT you must account for the impulses of coupled joints and contact impulses
        // by looping over the constraint of adjacent bodies.
        //
        for (int i = 0; i < numJoints; ++i)
        {

        }

        // TODO 
        // For each contact, compute an updated value of contacts[i]->lambda
        // using matrix-free pseudo-code provided in the notes.
        //
        // IMPORTANT you must account for the impulses of coupled joints and contact impulses
        // by looping over the constraint of adjacent bodies.
        //
        for(int i = 0; i < numContacts; ++i)
        {

        }
    }
}


