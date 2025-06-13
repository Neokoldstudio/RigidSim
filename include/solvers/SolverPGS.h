#pragma once

#include "solvers/Solver.h"
#include <vector>

// Projected Gauss-Seidel (PGS) Boxed constraint solver.
// This implementation uses a matrix-free approach, where only
// the non-zero blocks of the lead matrix are assembled and stored.
//
class SolverPGS : public Solver
{
public:

    SolverPGS(RigidBodySystem* _rigidBodySystem);

    // Implement PGS method that solves for the constraint impulses in @a m_rigidBodySystem.
    // 
    virtual void solve(float h) override;

private:

    std::vector<Eigen::MatrixXf> Ajoint;
    std::vector<Eigen::Matrix3f> Acontact;
    std::vector<Eigen::Vector3f> bcontact;
    std::vector<Eigen::VectorXf> bjoint;
};
