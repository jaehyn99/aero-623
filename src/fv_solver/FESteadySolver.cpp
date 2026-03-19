#include "FESteadySolver.h"
#include "CurvedElement.h"
#include "LocalTimeStepper.h"
#include "Residual.h"
#include "StateMesh.h"
#include "TimeIntegrator.h"
#include "TimeStepper.h"
#include "TriangularMesh.h"
#include <iostream>

FESteadySolver::FESteadySolver(std::shared_ptr<Residual> residual, std::shared_ptr<TimeIntegrator> integrator,
                               std::shared_ptr<TimeStepper> stepper, double tol):
    Solver(residual, integrator, stepper),
    _tol(tol)
{}

void FESteadySolver::solve(StateMesh& u) const{
    auto func = [this, &u](double t, Eigen::MatrixXd& x)
    {
        u.matrix() = std::move(x);
        auto mesh = u.mesh();
        const ReferenceElement& ref = mesh->reference();
        int Np = u.Np();

        Eigen::MatrixXd R = _residual->computeResidual(u);
        // std::cout << R.middleCols(0, Np) << std::endl << std::endl;
        // std::cout << mesh->elem(0).detJacobian() << std::endl;

        for (int k = 0; k < u.cellCount(); k++){
            Eigen::MatrixXd Rk = R.middleCols(k*Np, Np).transpose(); // create a temporary to avoid alias issues;
            if (!mesh->elem(k).isCurvedElement()){
                // On a linear element
                double detJ = mesh->elem(k).detJacobian();
                R.middleCols(k*Np, Np).transpose() = ref.MLLT().solve(Rk) / detJ;
            } else{
                // On a curved element
                CurvedElement& cElem = dynamic_cast<CurvedElement&>(mesh->elem(k));
                Eigen::LLT<Eigen::MatrixXd> MLLT = cElem.MLLT();
                R.middleCols(k*Np, Np).transpose() = ref.MLLT().solve(Rk);
            }
        }
        // std::cout << R.middleCols(0, Np) << std::endl << std::endl;
        return R;
    };

    double norm = func(0, u.matrix()).lpNorm<1>();
    _l1norm.push_back(norm); // the first L1-norm
    

    bool isConverged = false;
    int iter = 0;
    while (!isConverged){
        Eigen::ArrayXd dt = _stepper->dt(u);
        _integrator->integrate(func, u.matrix(), 0, dt);
        std::cout << u.cell(0) << std::endl << std::endl;
        //std::cout << u.cell(1) << std::endl << std::endl;
        assert(false);
        //u.computeGradient();
        norm = func(0, u.matrix()).lpNorm<1>();
        _l1norm.push_back(norm);
        isConverged = norm/_l1norm.front() <= _tol;
        iter++;
    }
    _result.emplace_back(u.matrix());
}