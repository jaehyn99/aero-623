#ifndef TRIDIAG_SOLVER_H
#define TRIDIAG_SOLVER_H

#include <Eigen/Dense>

namespace mesh {

// Solves tridiagonal system:
// A(i)*x(i) + C(i)*x(i+1) + B(i-1)*x(i-1) = D(i)
// with sizes:
// A: n, B: n-1, C: n-1, D: n
Eigen::VectorXd solveTridiag(const Eigen::VectorXd& A,
                            const Eigen::VectorXd& B,
                            const Eigen::VectorXd& C,
                            const Eigen::VectorXd& D);

} // namespace mesh

#endif
