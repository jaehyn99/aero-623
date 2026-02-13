#pragma once

#include <Eigen/Dense>

namespace solver {

using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using State = Eigen::Vector4d;
using Grad  = Eigen::Matrix<double,4,2>;

class SecondOrderEuler {
public:

    StateMatrix
    computeResidual(const Eigen::MatrixXi& I2E,
                    const Eigen::MatrixXi& B2E,
                    const Eigen::MatrixXd& In,
                    const Eigen::MatrixXd& Bn,
                    const Eigen::VectorXd& Area,
                    const StateMatrix& U) const;

private:

    Eigen::Vector4d computeNumericalFlux(const State& UL, const State& UR, const Eigen::Vector2d& normal) const;
    };
}; // namespace solver
