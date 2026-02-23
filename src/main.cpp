// src/main.cpp
#include <iostream>
#include <fstream>

#include "Constants.h"
#include "FVAdvectionFirstOrder.h"
#include "FVAdvectionSecondOrder.h"
#include "FVSteadySolver.h"
#include "FVUnsteadySolver.h"
#include "HLLEFlux.h"
#include "HybridWalkPNGrad.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "InviscidWallBC.h"
#include "LocalTimeStepper.h"
#include "OutletBC.h"
#include "RoeFlux.h"
#include "SSP_RK2.h"
#include "SSP_RK3.h"
#include "StateMesh.h"
#include "TriangularMesh.h"


int main() {
    std::shared_ptr<TriangularMesh> mesh;
    std::string meshName;
    do{
        std::cout << "Enter mesh name (\"coarse\", \"fine\", \"finer\", or \"finest\"): ";
        std::cin >> meshName;
        std::transform(meshName.begin(), meshName.end(), meshName.begin(), [](unsigned char c){ return std::tolower(c); });
    } while (meshName != "coarse" && meshName != "fine" && meshName != "finer" && meshName != "finest");
    if (meshName == "coarse") mesh = std::make_shared<TriangularMesh>("projects/Project-2/mesh_refined_2394.gri");
    else if (meshName == "fine") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined1.gri");
    else if (meshName == "finer") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined2.gri");
    else mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined3.gri");

    // Inlet conditions
    double gamma = 1.4;
    double rho0 = 1;
    double a0 = 1;
    double p0 = rho0*a0*a0/gamma;
    double alpha = 50*mconst::pi/180;
    double pout = 0.7*p0;
    double M = 0.1;

    // Boundary conditions
    std::shared_ptr<InletOutletBC> inlet = std::make_shared<InletOutletBC>(rho0, a0, alpha, pout, gamma);
    std::shared_ptr<BoundaryCondition> wall = std::make_shared<InviscidWallBC>(gamma);
    std::shared_ptr<BoundaryCondition> outlet = std::make_shared<OutletBC>(pout, gamma);
    std::vector<std::shared_ptr<BoundaryCondition>> bc{wall, inlet, wall, outlet};

    // Initialize the state mesh
    StateMesh states(mesh, bc);
    states.state(0).fill(rho0);
    states.state(1).fill(rho0*M*a0*std::cos(alpha));
    states.state(2).fill(rho0*M*a0*std::sin(alpha));
    states.state(3).fill(p0/(gamma-1) + 0.5*rho0*M*M*a0*a0);

    // Solver
    std::shared_ptr<FVFlux> flux;
    std::string fluxName;
    do{
        std::cout << "Enter flux name (\"roe\" or \"hlle\"): ";
        std::cin >> fluxName;
        std::transform(fluxName.begin(), fluxName.end(), fluxName.begin(), [](unsigned char c){ return std::tolower(c); });
    } while (fluxName != "roe" && fluxName != "hlle");
    if (fluxName == "roe") flux = std::make_shared<RoeFlux>(gamma);
    else flux = std::make_shared<HLLEFlux>(gamma);

    std::shared_ptr<FVResidual> residual;
    int FVOrder;
    do{
        std::cout << "Enter finite-volume order of accuracy (1 or 2): ";
        std::cin >> FVOrder;
    } while (FVOrder != 1 && FVOrder != 2);
    if (FVOrder == 1) residual = std::make_shared<FVAdvectionFirstOrder>(flux);
    else{
        states.setGradientMethod(std::make_shared<HybridWalkPNGrad>());
        int useLimiter;
        do{
            std::cout << "Will a BJ limiter be used (0 = No, 1 = Yes)";
            std::cin >> useLimiter;
        } while (useLimiter != 0 && useLimiter != 1);
        residual = std::make_shared<FVAdvectionSecondOrder>(flux, useLimiter==1);
    }
    
    std::shared_ptr<TimeIntegrator> integrator;
    int timeOrder;
    do{
        std::cout << "Enter time integration order of accuracy (2 or 3): ";
        std::cin >> timeOrder;    
    } while (timeOrder != 2 && timeOrder != 3);
    if (timeOrder == 2) integrator = std::make_shared<SSP_RK2>();
    else integrator = std::make_shared<SSP_RK3>();

    std::shared_ptr<TimeStepper> stepper = std::make_shared<LocalTimeStepper>(1, gamma, flux);
    std::unique_ptr<FVSolver> solver;
    int steadyState;
    int saveEveryNIterations, maxIterations;
    do{
        std::cout << "Enter solver mode (0 = unsteady, 1 = steady): ";
        std::cin >> steadyState;
    } while (steadyState != 0 && steadyState != 1);
    if (steadyState == 1) solver = std::make_unique<FVSteadySolver>(residual, integrator, stepper);
    else{
        inlet->setTransient(true);
        do{
            std::cout << "Enter the frequency (after every how many iterations) that data are saved: ";
            std::cin >> saveEveryNIterations;
        } while (saveEveryNIterations < 1);
        do{
            std::cout << "Enter the maximum number of iterations: ";
            std::cin >> maxIterations;
        } while (maxIterations < 1);
        solver = std::make_unique<FVUnSteadySolver>(residual, integrator, stepper, saveEveryNIterations, maxIterations);
    }

    try{
        if (FVOrder == 2 || steadyState == 0){
            // Creates a helper solver for second-order or unsteady simulations first
            std::shared_ptr<FVResidual> helperResidual = std::make_shared<FVAdvectionFirstOrder>(flux);
            FVSteadySolver helperSolver(helperResidual, integrator, stepper);
            helperSolver.solve(states);
        }

        solver->solve(states);
        std::vector<Eigen::MatrixXd> results = solver->getResult(); // size = 1 if steady, more than 1 if unsteady
        std::vector<double> l1norm = solver->getNorm();

        std::vector<std::size_t> iter;
        if (steadyState == 1) iter = {0};
        else{
            iter.resize(results.size());
            for (std::size_t i = 0; i < results.size(); i++) iter[i] = i;
        }

        std::ofstream file;
        std::string resultFilePath = "projects/Project-2/results/";
        resultFilePath += meshName + "_mesh_";
        resultFilePath += (steadyState == 0) ? "unsteady_" : "steady_";
        resultFilePath += (FVOrder == 1) ? "firstorder_" : "secondorder_";
        resultFilePath += (timeOrder == 2) ? "RK2_" : "RK3_";
        resultFilePath += fluxName;
        
        file.open(resultFilePath + "_norm.txt");
        for (auto norm: l1norm) file << norm << "\n";
        file.close();

        for (auto it: iter){
            if (steadyState == 0){
                std::string resultFilePathAtIter = resultFilePath + "_t_" + std::to_string(it*0.045*saveEveryNIterations);
                file.open(resultFilePathAtIter + ".txt");
            } else file.open(resultFilePath + ".txt");
            for (Eigen::Index e = 0; e < states.cellCount(); e++) file << results[it].col(e).transpose() << "\n";
            file.close();
        }
       
    } catch (std::runtime_error& ex){
        std::ofstream f("debug.txt");
    }

    return 0;
}