#pragma once

#include <Eigen/Dense>
#include "numericalFlux.hpp"
#include "boundaryState.hpp"
#include "boundaryFlux.hpp"
#include "mesh/TriangularMesh.h"

namespace solver {

class SecondOrderEuler {

private:
    const numericalFlux& numFlux_; // This can be Roe, HLLE, etc.
    const boundaryFlux& inletFlux_;
    const boundaryFlux& outletFlux_;
    const boundaryFlux& wallFlux_;
    const boundaryState& inletState_;
    const boundaryState& outletState_;
    const boundaryState& wallState_;
    double gamma_;

    using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;

    // Helper functions
    void computeGradient(TriangularMesh& mesh,
                         const Eigen::MatrixXi& I2E,
                         const Eigen::MatrixXi& B2E,
                         const Eigen::MatrixXd& In,
                         const Eigen::MatrixXd& Bn,
                         const Eigen::VectorXd& Area,
                         const StateMatrix& U,
                         StateMatrix& gradX,
                         StateMatrix& gradY) const;
    
    StateMatrix computeResidualFromGradient(TriangularMesh& mesh,
                                           const Eigen::MatrixXi& I2E,
                                           const Eigen::MatrixXi& B2E,
                                           const Eigen::MatrixXd& In,
                                           const Eigen::MatrixXd& Bn,
                                           const Eigen::VectorXd& Area,
                                           const StateMatrix& U,
                                           const StateMatrix& gradX,
                                           const StateMatrix& gradY) const;
    

public:
    SecondOrderEuler(const numericalFlux& numFl,
                      const boundaryFlux& inletFl,
                      const boundaryFlux& outletFl,
                      const boundaryFlux& wallFl,
                      const boundaryState& inletSt,
                      const boundaryState& outletSt,
                      const boundaryState& wallSt,
                      double gamma)
        : numFlux_(numFl), inletFlux_(inletFl), outletFlux_(outletFl), wallFlux_(wallFl),
          inletState_(inletSt), outletState_(outletSt), wallState_(wallSt),
          gamma_(gamma) {}

    StateMatrix computeResidual(TriangularMesh& mesh,
                                const Eigen::MatrixXi& I2E,
                                const Eigen::MatrixXi& B2E,
                                const Eigen::MatrixXd& In,
                                const Eigen::MatrixXd& Bn,
                                const Eigen::VectorXd& Area,
                                const StateMatrix& U) const;
};

} // namespace solver

