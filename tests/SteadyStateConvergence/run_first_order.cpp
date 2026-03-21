#include <iostream>
#include <Eigen/Dense>

#include "solver/FirstOrderEuler_alt.h"
#include "mesh/TriangularMesh.h"
#include "solver/hlleFlux.hpp"
#include "solver/inletFlux.hpp"
#include "solver/outletFlux.hpp"
#include "solver/wallFlux.hpp"
#include <fstream>
#include <vector>

// --------------------------------------------------
// Utility: load matrix from file
// --------------------------------------------------
Eigen::MatrixXd loadMatrix(const std::string& filename, int cols)
{
    std::ifstream file(filename);
    std::vector<double> values;
    double val;

    while (file >> val)
        values.push_back(val);

    int rows = values.size() / cols;
    Eigen::MatrixXd mat(rows, cols);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            mat(i,j) = values[i*cols + j];

    return mat;
}

Eigen::VectorXd loadVector(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<double> values;
    double val;

    while (file >> val)
        values.push_back(val);

    Eigen::VectorXd vec(values.size());
    for (size_t i = 0; i < values.size(); ++i)
        vec(i) = values[i];

    return vec;
}

int main()
{
    using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;

    // ------------------ Load Mesh ------------------
    TriangularMesh mesh("projects/Project-1/mesh_refined_2394.gri");

    Eigen::MatrixXi I2E = loadMatrix("projects/Project-1/mesh_refined_2394I2E.txt", 4).cast<int>();
    Eigen::MatrixXi B2E = loadMatrix("projects/Project-1/mesh_refined_2394B2E.txt", 3).cast<int>();
    Eigen::MatrixXd In  = loadMatrix("projects/Project-1/mesh_refined_2394In.txt", 2);
    Eigen::MatrixXd Bn  = loadMatrix("projects/Project-1/mesh_refined_2394Bn.txt", 2);
    Eigen::VectorXd Area = loadVector("projects/Project-1/mesh_refined_2394Area.txt");

    int nElem = Area.size();

    // ------------------ Initial Condition ------------------
    StateMatrix U(4, nElem);

    // Example: uniform freestream
    double rho = 0.791578;
    double u   = 0.402181;
    double v   = 0.47930;
    //double p   = 1.0;
    double gamma = 1.4;
    double R = 1.0;
    double rho0 = 1/R;
    double a0 = std::sqrt(gamma*R);

    double rhoE = 2.049597;

    for (int i = 0; i < nElem; ++i)
    {
        U(0,i) = rho;
        U(1,i) = rho*u;
        U(2,i) = rho*v;
        U(3,i) = rhoE;
    }

    // ------------------ Flux Objects ------------------
    HLLEFlux numFlux;
    inletFlux  inlet(rho0, a0, 0.8727, 0.0, false);
    outletFlux outlet(1.0);
    wallFlux   wall;

    solver::FirstOrderEuler_alt solver(
        numFlux,
        inlet,
        outlet,
        wall,
        gamma
    );

    // ------------------ March to Steady State ------------------
    solver.marchToSteadyState(mesh, I2E, B2E, In, Bn, Area, U, 0.1, 5000, 1e-6);

    std::cout << "Simulation finished.\n";

    return 0;
}
