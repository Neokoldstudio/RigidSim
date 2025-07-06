#include "rigidbody/RigidBodySystem.h"

#include "collision/CollisionDetect.h"
#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "solvers/SolverPGS.h"
#include <cstdio>

namespace
{
}

RigidBodySystem::RigidBodySystem() :
    solverIter(10), 
    m_preStepFunc(nullptr), 
    m_resetFunc(nullptr), 
    m_collisionsEnabled(true)
{
    m_collisionDetect = std::make_unique<CollisionDetect>(this);
    m_solver = new SolverPGS(this);
}

RigidBodySystem::~RigidBodySystem()
{
    clear();
}

void RigidBodySystem::addBody(RigidBody *_b)
{
    m_bodies.push_back(_b);
}

void RigidBodySystem::addJoint(Joint* _j)
{
    m_joints.push_back(_j);

    // Track list of joints for each body.
    // This makes it easy to find coupled constraints and contacts during the solve.
    if (_j->body0) _j->body0->joints.push_back(_j);
    if (_j->body1) _j->body1->joints.push_back(_j);
}

void RigidBodySystem::step(float dt)
{
    // Initialize the system.
    // Apply gravitional forces and reset angular forces.
    // Cleanup contacts from the previous time step.
    //
    for(auto b : m_bodies)
    {
        b->f = b->mass * Eigen::Vector3f(0.f, -9.81f, 0.f);
        b->tau.setZero();
        b->fc.setZero();
        b->tauc.setZero();
        b->contacts.clear();
    }

    // Standard simulation pipeline.
    computeInertias();

    if( m_preStepFunc )
    {
        m_preStepFunc(m_bodies);
    }

    m_collisionDetect->clear();

    if (m_collisionsEnabled)
    {
        m_collisionDetect->detectCollisions();
        m_collisionDetect->computeContactJacobians();
    }

    for (auto j : m_joints)
    {
        j->computeJacobian();
    }

    for(auto b : m_bodies)
    {
        b->fc.setZero();
        b->tauc.setZero();
    }

    calcConstraintForces(dt);

    for(RigidBody* b : m_bodies)
    {
        // TODO Objective #1
        // Numerical integration of rigid bodies.
        // For each rigid body in @m_bodies, perform numerical integration to:
        //   1. First update the velocities of each rigid body in @a m_bodies, 
        //   2. Followed by updates of the positions and orientations.
        //
        if (!b->fixed) {
            b->f += b->fc; // Add contact forces to the total force.
            b->tau += b->tauc; // Add contact torques to the total torque.
          
            Eigen::Matrix3f R = b->q.normalized().toRotationMatrix();
            b->Iinv = R * b->IbodyInv * R.transpose();

            // Update the linear and angular velocities.
            b->xdot += dt * (1 / b->mass) * b->f;
            b->omega += dt * b->Iinv * (b->tau - b->omega.cross(b->I*b->omega));

            // Update the position and orientation.
            b->x += dt * b->xdot;
            Eigen::Quaternionf dq;
            Eigen::Vector3f half_dt_omega = 0.5f * dt * b->omega;
            dq.w() = 1.0f;
            dq.vec() = half_dt_omega;
            b->q = (b->q * dq).normalized();
        }
    }
}

void RigidBodySystem::clear()
{
    if( m_resetFunc )
    {
        m_resetFunc();
    }

    m_collisionDetect->clear();

    // Cleanup all the joints.
    for (auto j : m_joints)
    {
        delete j;
    }
    m_joints.clear();

    // Finally, cleanup rigid bodies.
    for(auto b : m_bodies)
    {
        delete b;
    }
    m_bodies.clear();

}

void RigidBodySystem::computeInertias()
{
    for(RigidBody* b : m_bodies)
    {
        b->updateInertiaMatrix();
    }
}

// Accessors for the contact constraint array.
const std::vector<Contact>& RigidBodySystem::getContacts() const
{
    return m_collisionDetect->getContacts();
}

std::vector<Contact>& RigidBodySystem::getContacts()
{
    return m_collisionDetect->getContacts();
}

const std::vector<Joint*>& RigidBodySystem::getJoints() const
{
    return m_joints;
}

std::vector<Joint*>& RigidBodySystem::getJoints()
{
    return m_joints;
}

void RigidBodySystem::calcConstraintForces(float dt)
{
    // Solve for the constraint forces lambda
    //
    m_solver->setMaxIter(solverIter);
    m_solver->solve(dt);

    for (const auto j : m_joints)
    {
        const Eigen::Vector6f f0 = j->J0.transpose() * j->lambda / dt;
        const Eigen::Vector6f f1 = j->J1.transpose() * j->lambda / dt;

        j->body0->fc += f0.head<3>();
        j->body0->tauc += f0.tail<3>();
        j->body1->fc += f1.head<3>();
        j->body1->tauc += f1.tail<3>();
    }

    // Apply the constraint impulses as forces and torques acting on each body.
    // Essentially, we compute contact/joint forces by computing :
    //
    //       f = J^T * lambda / dt
    //
    // for each contact constraint.
    //
    auto contacts = m_collisionDetect->getContacts();
    for(const auto c : contacts)
    {
        // Convert impulses in c->lambda to forces.
        //
        const Eigen::Vector6f f0 = c.J0.transpose() * c.lambda / dt;
        const Eigen::Vector6f f1 = c.J1.transpose() * c.lambda / dt;

        // Apply forces to bodies that aren't fixed.
        if (!c.body0->fixed)
        {
            c.body0->fc += f0.head<3>();
            c.body0->tauc += f0.tail<3>();
        }
        if (!c.body1->fixed)
        {
            c.body1->fc += f1.head<3>();
            c.body1->tauc += f1.tail<3>();
        }
    }
}
