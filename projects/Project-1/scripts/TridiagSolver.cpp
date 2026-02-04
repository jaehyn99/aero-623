#include "mesh/TridiagSolver.h"

#include <stdexcept>

namespace mesh {

Eigen::VectorXd solveTridiag(const Eigen::VectorXd& A,
                            const Eigen::VectorXd& B,
                            const Eigen::VectorXd& C,
                            const Eigen::VectorXd& D)
{
    const Eigen::Index n = A.size();
    if (D.size() != n) throw std::invalid_argument("solveTridiag: D size mismatch");
    if (B.size() != n-1) throw std::invalid_argument("solveTridiag: B size mismatch");
    if (C.size() != n-1) throw std::invalid_argument("solveTridiag: C size mismatch");
    if (n < 2) throw std::invalid_argument("solveTridiag: n must be >= 2");

    Eigen::VectorXd CP = Eigen::VectorXd::Zero(n-1);
    Eigen::VectorXd DP = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd x  = Eigen::VectorXd::Zero(n);

    CP(0) = C(0) / A(0);
    DP(0) = D(0) / A(0);

    for (Eigen::Index i = 1; i < n-1; ++i) {
        const double denom = A(i) - B(i-1) * CP(i-1);
        CP(i) = C(i) / denom;
        DP(i) = (D(i) - B(i-1) * DP(i-1)) / denom;
    }

    // last row
    {
        const double denom = A(n-1) - B(n-2) * CP(n-2);
        DP(n-1) = (D(n-1) - B(n-2) * DP(n-2)) / denom;
    }

    // back substitution
    x(n-1) = DP(n-1);
    for (Eigen::Index k = n-2; k >= 0; --k) {
        x(k) = DP(k) - CP(k) * x(k+1);
        if (k == 0) break; // avoid underflow of signed index
    }

    return x;
}

} // namespace mesh
