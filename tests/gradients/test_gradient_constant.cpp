// #include <iostream>
// #include <cassert>
// #include <fstream>
// #include <Eigen/Dense>

// #include "solver/SecondorderEuler.h"
// #include "mesh/TriangularMesh.h"
// #include "ConstantBoundaryState.h"

// // --------------------------------------------------
// // Utility: load matrix from file
// // --------------------------------------------------
// Eigen::MatrixXd loadMatrix(const std::string& filename, int cols)
// {
//     std::ifstream file(filename);
//     std::vector<double> values;
//     double val;

//     while (file >> val)
//         values.push_back(val);

//     int rows = values.size() / cols;
//     Eigen::MatrixXd mat(rows, cols);

//     for (int i = 0; i < rows; ++i)
//         for (int j = 0; j < cols; ++j)
//             mat(i,j) = values[i*cols + j];

//     return mat;
// }

// Eigen::VectorXd loadVector(const std::string& filename)
// {
//     std::ifstream file(filename);
//     std::vector<double> values;
//     double val;

//     while (file >> val)
//         values.push_back(val);

//     Eigen::VectorXd vec(values.size());
//     for (size_t i = 0; i < values.size(); ++i)
//         vec(i) = values[i];

//     return vec;
// }

// // --------------------------------------------------
// // Dummy fluxes (not used for gradient test)
// // --------------------------------------------------
// class DummyFlux : public numericalFlux {
// public:
//     Eigen::Vector4d operator()(
//         const Eigen::Vector4d&,
//         const Eigen::Vector4d&,
//         double,
//         const Eigen::Vector2d&) const override
//     {
//         return Eigen::Vector4d::Zero();
//     }
// };

// class DummyBoundaryFlux : public boundaryFlux {
// public:
//     Eigen::Vector4d operator()(
//         const Eigen::Vector4d&,
//         double,
//         const Eigen::Vector2d&) const override
//     {
//         return Eigen::Vector4d::Zero();
//     }
// };

// // --------------------------------------------------
// // MAIN TEST
// // --------------------------------------------------
// int main()
// {
//     using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;

//     // --------------------------------------------------
//     // Load mesh
//     // --------------------------------------------------
//     TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");

//     // Generate connectivity files
//     mesh.writeGri("projects/Project-1/mesh_coarse");

//     // --------------------------------------------------
//     // Load connectivity data
//     // --------------------------------------------------
//     Eigen::MatrixXi I2E = loadMatrix("projects/Project-1/mesh_coarseI2E.txt", 4).cast<int>();
//     Eigen::MatrixXi B2E = loadMatrix("projects/Project-1/mesh_coarseB2E.txt", 3).cast<int>();
//     Eigen::MatrixXd In  = loadMatrix("projects/Project-1/mesh_coarseIn.txt", 2);
//     Eigen::MatrixXd Bn  = loadMatrix("projects/Project-1/mesh_coarseBn.txt", 2);
//     Eigen::VectorXd Area = loadVector("projects/Project-1/mesh_coarseArea.txt");

//     // debug line
//     std::cout << I2E.minCoeff() << " " << I2E.maxCoeff() << std::endl;
//     std::cout << B2E.minCoeff() << " " << B2E.maxCoeff() << std::endl;


//     // --------------------------------------------------
//     // Create solver with CONSTANT boundary state
//     // --------------------------------------------------
//     DummyFlux numFlux;
//     DummyBoundaryFlux inletFlux;
//     DummyBoundaryFlux outletFlux;
//     DummyBoundaryFlux wallFlux;

//     ConstantBoundaryState inletState;
//     ConstantBoundaryState outletState;
//     ConstantBoundaryState wallState;

//     solver::SecondOrderEuler solver(numFlux, inletFlux, outletFlux, wallFlux, inletState, outletState, wallState, 1.4);

//     // --------------------------------------------------
//     // Constant state test
//     // --------------------------------------------------
//     int nElem = Area.size();

//     StateMatrix U = StateMatrix::Ones(4, nElem);
//     StateMatrix gradX = StateMatrix::Zero(4, nElem);
//     StateMatrix gradY = StateMatrix::Zero(4, nElem);

//     solver.test_ComputeGradient(mesh, I2E, B2E, In, Bn, Area, U, gradX, gradY);

//     double tol = 1e-12;

//     for (int e = 0; e < nElem; ++e)
//         for (int v = 0; v < 4; ++v)
//             assert(std::abs(gradX(v,e)) < tol &&
//                    std::abs(gradY(v,e)) < tol);

//     std::cout << "Gradient constant-state test PASSED\n";

//     return 0;
// }


#include <iostream>
#include <cassert>
#include <fstream>
#include <set>
#include <Eigen/Dense>

#include "solver/SecondorderEuler.h"
#include "mesh/TriangularMesh.h"
#include "ConstantBoundaryState.h"

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

// --------------------------------------------------
// Write VTK file for ParaView
// --------------------------------------------------
void writeVTK(const TriangularMesh& mesh,
              const Eigen::MatrixXd& gradX,
              const Eigen::MatrixXd& gradY,
              const std::string& filename)
{
    std::ofstream out(filename);

    int nNodes = mesh.numNodes();
    int nElems = mesh.numElems();

    out << "# vtk DataFile Version 3.0\n";
    out << "Gradient Output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // ------------------ POINTS ------------------
    out << "POINTS " << nNodes << " double\n";
    for (int i = 0; i < nNodes; ++i)
    {
        auto p = mesh.node(i);
        out << p(0) << " " << p(1) << " 0.0\n";
    }

    // ------------------ CELLS ------------------
    out << "CELLS " << nElems << " " << 4*nElems << "\n";
    for (int e = 0; e < nElems; ++e)
    {
        auto elem = mesh.elem(e);
        out << "3 "
            << elem._pointID[0] << " "
            << elem._pointID[1] << " "
            << elem._pointID[2] << "\n";
    }

    // ------------------ CELL TYPES ------------------
    out << "CELL_TYPES " << nElems << "\n";
    for (int e = 0; e < nElems; ++e)
        out << "5\n";  // VTK_TRIANGLE

    // ------------------ CELL DATA ------------------
    out << "CELL_DATA " << nElems << "\n";

    for (int v = 0; v < 4; ++v)
    {
        out << "SCALARS gradX_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e)
            out << gradX(v,e) << "\n";

        out << "SCALARS gradY_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e)
            out << gradY(v,e) << "\n";
    }

    out.close();
}

// --------------------------------------------------
// Dummy fluxes (not used for gradient test)
// --------------------------------------------------
class DummyFlux : public numericalFlux {
public:
    Eigen::Vector4d operator()(
        const Eigen::Vector4d&,
        const Eigen::Vector4d&,
        double,
        const Eigen::Vector2d&) const override
    {
        return Eigen::Vector4d::Zero();
    }
};

class DummyBoundaryFlux : public boundaryFlux {
public:
    Eigen::Vector4d operator()(
        const Eigen::Vector4d&,
        double,
        const Eigen::Vector2d&) const override
    {
        return Eigen::Vector4d::Zero();
    }
};

// --------------------------------------------------
// MAIN TEST
// --------------------------------------------------
int main()
{
    using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;

    // ------------------ Load Mesh ------------------
    TriangularMesh mesh("projects/Project-1/mesh_refined_2394.gri");
    mesh.writeGri("projects/Project-1/mesh_refined_2394");

    // ------------------ Load Connectivity ------------------
    Eigen::MatrixXi I2E = loadMatrix("projects/Project-1/mesh_refined_2394I2E.txt", 4).cast<int>();
    Eigen::MatrixXi B2E = loadMatrix("projects/Project-1/mesh_refined_2394B2E.txt", 3).cast<int>();
    Eigen::MatrixXd In  = loadMatrix("projects/Project-1/mesh_refined_2394In.txt", 2);
    Eigen::MatrixXd Bn  = loadMatrix("projects/Project-1/mesh_refined_2394Bn.txt", 2);
    Eigen::VectorXd Area = loadVector("projects/Project-1/mesh_refined_2394Area.txt");

    std::cout << "I2E rows: " << I2E.rows() << "\n";
    std::cout << "In rows:  " << In.rows()  << "\n";
    std::cout << "B2E rows: " << B2E.rows() << "\n";
    std::cout << "Bn rows:  " << Bn.rows()  << "\n";

    // ------------------ Solver Setup ------------------
    DummyFlux numFlux;
    DummyBoundaryFlux inletFlux;
    DummyBoundaryFlux outletFlux;
    DummyBoundaryFlux wallFlux;

    ConstantBoundaryState inletState;
    ConstantBoundaryState outletState;
    ConstantBoundaryState wallState;

    solver::SecondOrderEuler solver(
        numFlux,
        inletFlux,
        outletFlux,
        wallFlux,
        inletState,
        outletState,
        wallState,
        1.4);

    // ------------------ Constant State Test ------------------
    int nElem = Area.size();

    StateMatrix U = StateMatrix::Ones(4, nElem);
    StateMatrix gradX = StateMatrix::Zero(4, nElem);
    StateMatrix gradY = StateMatrix::Zero(4, nElem);

    solver.test_ComputeGradient(mesh, I2E, B2E, In, Bn, Area, U, gradX, gradY);

    // ------------------ Write VTK ------------------
    writeVTK(mesh, gradX, gradY, "gradient_output.vtk");

    std::cout << "VTK file written: gradient_output.vtk\n";

    // double tol = 1e-12;

    // for (int e = 0; e < nElem; ++e)
    //     for (int v = 0; v < 4; ++v)
    //         assert(std::abs(gradX(v,e)) < tol &&
    //                std::abs(gradY(v,e)) < tol);

    // std::cout << "Gradient constant-state test PASSED\n";

    

    return 0;
}
