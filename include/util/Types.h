#pragma once

#include <Eigen/Dense>

// Define some important matrix, vector types.
// This file can also be used to override default scalar types
//   e.g. for autodiff computations or changing fp precision
//
namespace Eigen
{
	// Vector types	
	typedef Eigen::Matrix< float, 6, 1 > Vector6f;

	// Matrix types	
	typedef Eigen::Matrix< float, 6, 6 > Matrix6f;
}

typedef Eigen::Matrix<float, Eigen::Dynamic, 6, Eigen::RowMajor> JBlock;		// generic constraint Jacobian
typedef Eigen::Matrix<float, Eigen::Dynamic, 6, Eigen::ColMajor> JBlockTranspose;	// constraint Jacobian transpose

typedef Eigen::Matrix<float, 3, 6, Eigen::RowMajor> JCBlock;			// contact Jacobian
typedef Eigen::Matrix<float, 3, 6, Eigen::ColMajor> JCBlockTranspose;	// contact Jacobian transpose
