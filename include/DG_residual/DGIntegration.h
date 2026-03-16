#ifndef DG_INTEGRATION_H
#define DG_INTEGRATION_H

#include <Eigen/Dense>

class DGIntegration
{
private:

public:
    virtual ~DGIntegration() = default;
    virtual DGIntegration() = default;

    void DGIntegration::computeResidual(int p_,
                                        double gamma,
                                        const TriangularMesh& mesh,
                                        const std::vector<std::shared_ptr<BoundaryCondition>>& bc,
                                        const FVFlux& flux,
                                        const Eigen::MatrixXd& modes,
                                        Eigen::MatrixXd& residual);

    void getElemIntegral(int p_,
                        double gamma,
                        const TriangularMesh& mesh,
                        const Eigen::MatrixXd& modes,
                        Eigen::MatrixXd& residual);

    void getLineIntegral(int p_,
                         double gamma,
                         const TriangularMesh& mesh,
                         const std::vector<std::shared_ptr<BoundaryCondition>>& bc,
                         const FVFlux& flux,
                         const Eigen::MatrixXd& modes,
                         Eigen::MatrixXd& residual);


};

