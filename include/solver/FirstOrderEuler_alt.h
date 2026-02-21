#pragma once

#include <Eigen/Dense>
#include "numericalFlux.hpp"
#include "boundaryState.hpp"
#include "boundaryFlux.hpp"
#include "mesh/TriangularMesh.h"

using State       = Eigen::Vector4d;
using StateMatrix = Eigen::Matrix<double, 4, Eigen::Dynamic>;

struct ResidualResult {
    StateMatrix R;
    Eigen::VectorXd waveSpeed;
    Eigen::VectorXd perimeter;
};

namespace solver {

class FirstOrderEuler_alt {

private:
    const numericalFlux& numFlux_; // This can be Roe, HLLE, etc.
    const boundaryFlux& inletFlux_;
    const boundaryFlux& outletFlux_;
    const boundaryFlux& wallFlux_;
    double gamma_;

public:
    FirstOrderEuler_alt(const numericalFlux& numFl,
                         const boundaryFlux& inletFl,
                         const boundaryFlux& outletFl,
                         const boundaryFlux& wallFl,
                        double gamma)
        : numFlux_(numFl), inletFlux_(inletFl), outletFlux_(outletFl), wallFlux_(wallFl),
          gamma_(gamma) {}
    
    // ======================================================
    // Public Solver
    // ======================================================
    ResidualResult computeResidual(TriangularMesh& mesh,
                                   const Eigen::MatrixXi& I2E,
                                   const Eigen::MatrixXi& B2E,
                                   const Eigen::MatrixXd& In,
                                   const Eigen::MatrixXd& Bn,
                                   const Eigen::VectorXd& Area,
                                   const StateMatrix& U) const;
    
    void marchToSteadyState(TriangularMesh& mesh,
                            const Eigen::MatrixXi& I2E,
                            const Eigen::MatrixXi& B2E,
                            const Eigen::MatrixXd& In,
                            const Eigen::MatrixXd& Bn,
                            const Eigen::VectorXd& Area,
                            StateMatrix& U,
                            double CFL,
                            int maxIter,
                            double tol) const;

};

} // namespace solver

