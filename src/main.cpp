// src/main.cpp
#include <iostream>
#include <fstream>

#include "Constants.h"
#include "Element.h"
#include "Face.h"
#include "FEAdvection.h"
#include "FESteadySolver.h"
#include "GaussLegendre1D.h"
#include "GaussLegendre2D.h"
#include "HLLEFlux.h"
#include "InletBC.h"
#include "InletOutletBC.h"
#include "InviscidWallBC.h"
#include "LocalTimeStepper.h"
#include "OutletBC.h"
#include "RK4.h"
#include "RoeFlux.h"
#include "SSP_RK3.h"
#include "StateMesh.h"
#include "TriangularMesh.h"

#include "Lagrange2DBasisFunctions.h"

int main() {
    std::shared_ptr<TriangularMesh> mesh;
    std::string meshName;

    // std::cout << std::endl;
    // for (int k = 0; k < int(mesh->numElems()); k++){
    //     std::cout << "Element " << k << " has these points: " << std::endl;
    //     const Element& elem = mesh->elem(k);
    //     for (int i = 0; i < 3; i++){
    //         std::cout << "\tPoint " << elem.pointID(i) << " at [" << mesh->node(elem.pointID(i)).transpose() << "]." << std::endl; 
    //     }

    //     std::cout << "Element " << k << " has these edges: " << std::endl;
    //     for (int i = 0; i < 3; i++){
    //         const Face& face = mesh->face(elem.faceID(i));
    //         std::cout << "\tEdge " << elem.faceID(i) << " containing points " << face.pointID(0) << " and " << face.pointID(1) << ". ";
    //         if (face.isBoundaryFace()){
    //             std::cout << "It is on boundary " << face.title();
    //             std::cout << ". Normal vector = " << face.normal().transpose();
    //         } else if (face.isPeriodicFace()){
    //             std::cout << "It is periodic with element " << (face.elemID(0) == k ? face.elemID(1) : face.elemID(0));
    //             std::cout << ". Normal vector = " << face.normal().transpose() << " pointing from " << face.elemID(0) << " to " << face.elemID(1);
    //         } else{
    //             std::cout << "It is an internal edge and neighbors element " << (face.elemID(0) == k ? face.elemID(1) : face.elemID(0));
    //             std::cout << ". Normal vector = " << face.normal().transpose() << " pointing from " << face.elemID(0) << " to " << face.elemID(1); 
    //         }
    //         std::cout << "." << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    do{
        std::cout << "Enter mesh name (\"test\", \"coarse\", \"fine\", \"finer\", or \"finest\"): ";
        std::cin >> meshName;
        std::transform(meshName.begin(), meshName.end(), meshName.begin(), [](unsigned char c){ return std::tolower(c); });
    } while (meshName != "test" && meshName != "coarse" && meshName != "fine" && meshName != "finer" && meshName != "finest");

    std::size_t p = 0; // Lagrange order for solution approx
    std::size_t q = 1; // Lagrange order for geometry approx
    // do{
    //     std::cout << "Enter the Lagrange polynomial order for solution approximation (p = 0, 1, 2, or 3): ";
    //     std::cin >> p;
    // } while (p < 0 || p > 3);
    // do{
    //     std::cout << "Enter the Lagrange polynomial order for geometry approximation (q = 1 or 3): ";
    //     std::cin >> q;
    // } while (q != 1 && q != 3);    
    std::size_t r = 2*(p+q);

    if (meshName == "test") mesh = std::make_shared<TriangularMesh>("projects/Project-3/test2.gri", p, q, r, false);
    else if (meshName == "coarse") mesh = std::make_shared<TriangularMesh>("projects/Project-2/mesh_refined_2394.gri", p, q, r);
    else if (meshName == "fine") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined1.gri", p, q, r);
    else if (meshName == "finer") mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined2.gri", p, q, r);
    else mesh = std::make_shared<TriangularMesh>("projects/Project-2/meshGlobalRefined3.gri", p, q, r);

    // Inlet conditions
    double gamma = 1.4;
    double rho0 = 1;
    double a0 = 1;
    double p0 = rho0*a0*a0/gamma;
    double alpha = 50*mconst::pi/180;
    double pout = 0.7*p0;
    double M = 0.1;

    // Boundary conditions
    // std::shared_ptr<InletOutletBC> inlet = std::make_shared<InletOutletBC>(rho0, a0, alpha, pout, gamma);
    // std::shared_ptr<BoundaryCondition> wall = std::make_shared<InviscidWallBC>(gamma);
    // std::shared_ptr<BoundaryCondition> outlet = std::make_shared<OutletBC>(pout, gamma);
    // std::vector<std::shared_ptr<BoundaryCondition>> bc{wall, inlet, wall, outlet};
    std::shared_ptr<InletBC> inlet = std::make_shared<InletBC>(rho0, a0, alpha, gamma);
    std::vector<std::shared_ptr<BoundaryCondition>> bc{nullptr, nullptr, nullptr, nullptr};

    // Initialize the state mesh
    double rhoi = rho0;
    double rhoui = rho0*M*a0*std::cos(alpha);
    double rhovi = rho0*M*a0*std::sin(alpha);
    double rhoEi = p0/(gamma-1) + 0.5*rho0*M*M*a0*a0;

    StateMesh U(mesh, bc, 4, p);
    U.state(0).fill(rho0);
    U.state(1).fill(rho0*M*a0*std::cos(alpha));
    U.state(2).fill(rho0*M*a0*std::sin(alpha));
    U.state(3).fill(p0/(gamma-1) + 0.5*rho0*M*M*a0*a0);

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

    std::shared_ptr<Residual> residual = std::make_shared<FEAdvection>(flux);
    std::cout << residual->computeResidual(U) << std::endl;

    // // int FVOrder;
    // // do{
    // //     std::cout << "Enter finite-volume order of accuracy (1 or 2): ";
    // //     std::cin >> FVOrder;
    // // } while (FVOrder != 1 && FVOrder != 2);
    // // if (FVOrder == 1) residual = std::make_shared<FVAdvectionFirstOrder>(flux);
    // // else{
    // //     U.setGradientMethod(std::make_shared<HybridWalkPNGrad>());
    // //     int useLimiter;
    // //     do{
    // //         std::cout << "Will a BJ limiter be used (0 = No, 1 = Yes)";
    // //         std::cin >> useLimiter;
    // //     } while (useLimiter != 0 && useLimiter != 1);
    // //     residual = std::make_shared<FVAdvectionSecondOrder>(flux, useLimiter==1);
    // // }
    
    std::shared_ptr<TimeIntegrator> integrator;
    int timeOrder;
    do{
        std::cout << "Enter time integration order of accuracy (3 or 4): ";
        std::cin >> timeOrder;    
    } while (timeOrder != 3 && timeOrder != 4);
    if (timeOrder == 3) integrator = std::make_shared<SSP_RK3>();
    else integrator = std::make_shared<RK4>();

    std::shared_ptr<TimeStepper> stepper = std::make_shared<LocalTimeStepper>(1, gamma, flux);
    std::unique_ptr<Solver> solver = std::make_unique<FESteadySolver>(residual, integrator, stepper);
    int steadyState = 1;

    // int saveEveryNIterations, maxIterations;
    // do{
    //     std::cout << "Enter solver mode (0 = unsteady, 1 = steady): ";
    //     std::cin >> steadyState;
    // } while (steadyState != 0 && steadyState != 1);
    // if (steadyState == 1) solver = std::make_unique<FVSteadySolver>(residual, integrator, stepper);
    // else{
    //     inlet->setTransient(true);
    //     do{
    //         std::cout << "Enter the frequency (after every how many iterations) that data are saved: ";
    //         std::cin >> saveEveryNIterations;
    //     } while (saveEveryNIterations < 1);
    //     do{
    //         std::cout << "Enter the maximum number of iterations: ";
    //         std::cin >> maxIterations;
    //     } while (maxIterations < 1);
    //     solver = std::make_unique<FVUnSteadySolver>(residual, integrator, stepper, saveEveryNIterations, maxIterations);
    // }

    try{
        // if (FVOrder == 2 || steadyState == 0){
        //     // Creates a helper solver for second-order or unsteady simulations first
        //     std::shared_ptr<Residual> helperResidual = std::make_shared<FVAdvectionFirstOrder>(flux);
        //     FVSteadySolver helperSolver(helperResidual, integrator, stepper);
        //     helperSolver.solve(U);
        // }

        solver->solve(U);
        std::vector<Eigen::MatrixXd> results = solver->getResult(); // size = 1 if steady, more than 1 if unsteady
        for (auto& Ui: results){
            Ui.row(0).array() -= rhoi;
            Ui.row(1).array() -= rhoui;
            Ui.row(2).array() -= rhovi;
            Ui.row(3).array() -= rhoEi;
        }
        std::cout << results[0] << std::endl << std::endl;
        std::cout << results[1] << std::endl << std::endl;
        std::cout << results[2] << std::endl << std::endl;
        
        // std::vector<double> l1norm = solver->getNorm();

        // std::vector<std::size_t> iter;
        // if (steadyState == 1) iter = {0};
        // else{
        //     iter.resize(results.size());
        //     for (std::size_t i = 0; i < results.size(); i++) iter[i] = i;
        // }

        // std::ofstream file;
        // std::string resultFilePath = "projects/Project-3/results/";
        // resultFilePath += meshName + "_mesh_";
        // resultFilePath += (steadyState == 0) ? "unsteady_" : "steady_";
        // resultFilePath += "p" + std::to_string(p) + "_";
        // resultFilePath += "q" + std::to_string(q) + "_";
        // resultFilePath += (timeOrder == 3) ? "RK3_" : "RK4_";
        // resultFilePath += fluxName;
        
        // file.open(resultFilePath + "_norm.txt");
        // for (auto norm: l1norm) file << norm << "\n";
        // file.close();

        // for (auto it: iter){
        //     if (steadyState == 0){
        //         std::string resultFilePathAtIter = resultFilePath + "_t_" + std::to_string(it*0.045*saveEveryNIterations);
        //         file.open(resultFilePathAtIter + ".txt");
        //     } else file.open(resultFilePath + ".txt");
        //     for (Eigen::Index e = 0; e < U.cellCount(); e++) file << results[it].col(e).transpose() << "\n";
        //     file.close();
        // }
       
    } catch (std::runtime_error& ex){
        std::cerr << ex.what() << std::endl;
    }

    return 0;
}