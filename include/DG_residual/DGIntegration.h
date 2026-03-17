#ifndef DG_INTEGRATION_H
#define DG_INTEGRATION_H

#include "Eigen/Dense"

class TriangularMesh;
class FVFlux;
class BoundaryCondition;
class DGIntegration{
    public:
    DGIntegration(int p, int q);
    virtual ~DGIntegration() = default;

    void DGIntegration::computeResidual(double gamma,
                                        const TriangularMesh& mesh,
                                        const std::vector<std::shared_ptr<BoundaryCondition>>& bc,
                                        const FVFlux& flux,
                                        const Eigen::MatrixXd& modes,
                                        Eigen::MatrixXd& residual);

    // Compute area integral of grad(phi) dot F
    void getElemIntegral(double gamma,
                        const TriangularMesh& mesh,
                        const Eigen::MatrixXd& modes,
                        Eigen::MatrixXd& residual);

    // Compute line integral of phi F
    void getLineIntegral(double gamma,
                        const TriangularMesh& mesh,
                        const std::vector<std::shared_ptr<BoundaryCondition>>& bc,
                        const FVFlux& flux,
                        const Eigen::MatrixXd& modes,
                        Eigen::MatrixXd& residual);
        
    protected:
    int _p;
    int _q;
};

#endif

