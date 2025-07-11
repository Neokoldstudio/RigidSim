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
        Eigen::Vector3f r = b - A * x;
        Eigen::Vector3f dx = A.ldlt().solve(r);
        x += dx;

        x(0) = std::max(0.0f, x(0));

        float friction_limit = mu * x(0);
        Eigen::Vector2f t(x(1), x(2));
        float tNorm = t.norm();
        float epsilon = 1e-8f;
        if (tNorm > friction_limit) {// the friction impulse is too large, it needs to be projected back onto the friction cone
            float scale = friction_limit / (tNorm + epsilon);
            x(1) *= scale;
            x(2) *= scale;
        }
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
        // Solve the linear system A * dx = b - A * x for dx, then update x
        Eigen::VectorXf r = b - A * x;
        Eigen::VectorXf dx = A.ldlt().solve(r);
        x += dx;
    }

    void accumulateCoupledContactsAndJoints(Joint* j, RigidBody* body, Eigen::VectorXf& b)
    {
        // Accumulate contributions from contacts
        for (Contact* c : body->contacts)
        {
            if (c->body0 == body)
            {
                b += c->J0Minv * c->J0.transpose() * c->lambda;
            }
            else if (c->body1 == body)
            {
                b += c->J1Minv * c->J1.transpose() * c->lambda;
            }
        }
    
        // Accumulate contributions from other joints
        for (Joint* j2 : body->joints)
        {
            if (j2 == j) continue;
    
            if (j2->body0 == body)
            {
                b += j2->J0Minv * j2->J0.transpose() * j2->lambda;
            }
            else if (j2->body1 == body)
            {
                b += j2->J1Minv * j2->J1.transpose() * j2->lambda;
            }
        }
    }

    void accumulateCoupledContactsAndJoints(Contact* c, RigidBody* body, Eigen::Vector3f& b)
    {
        // Accumulate contributions from contacts
        for (Contact *c2 : body->contacts)
        {
            if (c2 == c) continue; // Skip the current contact
            
            if (c->body0 == body)
            {
                b += c2->J0Minv * c2->J0.transpose() * c2->lambda;
            }
            else if (c2->body1 == body)
            {
                b += c2->J1Minv * c2->J1.transpose() * c2->lambda;
            }
        }
    
        // Accumulate contributions from other joints
        for (Joint* j : body->joints)
        {
    
            if (j->body0 == body)
            {
                b += j->J0Minv * j->J0.transpose() * j->lambda;
            }
            else if (j->body1 == body)
            {
                b += j->J1Minv * j->J1.transpose() * j->lambda;
            }
        }
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
            Contact* c = &contacts[i];
            const float eps = 1e-6f;

            // Compute the diagonal term: Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            Eigen::Matrix3f A = eps * Eigen::Matrix3f::Identity();

            if (!c->body0->fixed)
            {
                A += c->J0Minv * c->J0.transpose();
            }
            if (!c->body1->fixed)
            {
                A += c->J1Minv * c->J1.transpose();
            }

            Acontact[i] = A;
        }
    }

    bcontact.resize(numContacts);
    for (int i = 0; i < numContacts; ++i)
    {
        Contact* c = &contacts[i];
        c->lambda.setOnes();

        // TODO Compute the right-hand side vector for contacts, e.g.
        //      b = -gamma*phi/h - J*vel - h*JMinvJT*force
        //
        float gamma = 0.7f;

        Eigen::Matrix<float, 6, 1> b0Force, b1Force;
        b0Force << c->body0->f, c->body0->tau;
        b1Force << c->body1->f, c->body1->tau;
    
        Eigen::Matrix<float, 6, 1> b0Vel, b1Vel;
        b0Vel << c->body0->xdot, c->body0->omega;
        b1Vel << c->body1->xdot, c->body1->omega;
    
        Eigen::Vector3f b = 
            -gamma * (1.0 / h)* c->phi
            - (c->J0 * b0Vel + c->J1 * b1Vel) 
            - h * (c->J0Minv * b0Force + c->J1Minv * b1Force);
    
        bcontact[i] = b;
    }

    bjoint.resize(numJoints);
    for (int i = 0; i < numJoints; ++i)
    {
        Joint* j = joints[i];

        // TODO Compute the right-hand side vector for joints, e.g.
        //      b = -gamma*phi/h - J*vel - h*JMinvJT*force
        //
        float gamma = 0.7f;
        Eigen::Matrix<float, 6, 1> b0Force, b1Force;
        b0Force << j->body0->f, j->body0->tau;
        b1Force << j->body1->f, j->body1->tau;

        Eigen::Matrix<float, 6, 1> b0Vel, b1Vel;
        b0Vel << j->body0->xdot, j->body0->omega;
        b1Vel << j->body1->xdot, j->body1->omega;

        Eigen::VectorXf b = -gamma * (1.0 / h) * j->phi -
                            (j->J0 * b0Vel + j->J1 * b1Vel) -
                            h * (j->J0Minv * b0Force + j->J1Minv * b1Force);
        
        bjoint[i] = b;
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
            Joint *j = joints[i];
            // Accumulate contributions from contacts and other joints
            auto b = bjoint[i];
            accumulateCoupledContactsAndJoints(j, j->body0, b);
            accumulateCoupledContactsAndJoints(j, j->body1, b);
            solveJoint(Ajoint[i], bjoint[i], j->lambda);
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
            Contact *c = &contacts[i];
            // Accumulate contributions from contacts and other joints
            Eigen::Vector3f b = bcontact[i];
            accumulateCoupledContactsAndJoints(c, c->body0, b);
            accumulateCoupledContactsAndJoints(c, c->body1, b);
        
            solveContact(Acontact[i], bcontact[i], c->lambda, c->mu);
        }
    }
}


