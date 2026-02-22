
// tests/test_gradient_constant.cpp
#include <iostream>
#include <cassert>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <cmath>
#include <vector>
#include <memory>
#include <algorithm>

#include "mesh/TriangularMesh.h"
#include "mesh/StateMesh.h"
#include "ConstantBC.h"                      
#include "fv_FVGradient/hybridWalkPNGrad.h"
#include "fv_FVGradient/WalkGrad.h"

using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using GradCell = Eigen::Matrix<double,4,2>;
using GradVec  = std::vector<GradCell, Eigen::aligned_allocator<GradCell>>;

// --------------------------------------------------
// Utility: load matrix from file
// --------------------------------------------------
Eigen::MatrixXd loadMatrix(const std::string& filename, int cols)
{
    std::ifstream file(filename);
    std::vector<double> values;
    double val;
    while (file >> val) values.push_back(val);

    int rows = int(values.size()) / cols;
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
    while (file >> val) values.push_back(val);

    Eigen::VectorXd vec(values.size());
    for (size_t i = 0; i < values.size(); ++i) vec(i) = values[i];
    return vec;
}

// --------------------------------------------------
// Write VTK file for ParaView
// --------------------------------------------------
void writeVTK(const TriangularMesh& mesh,
              const StateMatrix& U,
              const Eigen::MatrixXd& gradX,
              const Eigen::MatrixXd& gradY,
              const std::string& filename)
{
    std::ofstream out(filename);

    int nNodes = int(mesh.numNodes());
    int nElems = int(mesh.numElems());

    out << "# vtk DataFile Version 3.0\n";
    out << "State and Gradient Output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << nNodes << " double\n";
    for (int i = 0; i < nNodes; ++i)
    {
        auto p = mesh.node(i);
        out << p(0) << " " << p(1) << " 0.0\n";
    }

    out << "CELLS " << nElems << " " << 4*nElems << "\n";
    for (int e = 0; e < nElems; ++e)
    {
        auto elem = mesh.elem(e);
        out << "3 "
            << elem._pointID[0] << " "
            << elem._pointID[1] << " "
            << elem._pointID[2] << "\n";
    }

    out << "CELL_TYPES " << nElems << "\n";
    for (int e = 0; e < nElems; ++e) out << "5\n";

    out << "CELL_DATA " << nElems << "\n";

    for (int v = 0; v < 4; ++v)
    {
        out << "SCALARS U_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e) out << U(v,e) << "\n";
    }

    for (int v = 0; v < 4; ++v)
    {
        out << "SCALARS gradX_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e) out << gradX(v,e) << "\n";
    }

    for (int v = 0; v < 4; ++v)
    {
        out << "SCALARS gradY_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e) out << gradY(v,e) << "\n";
    }
}

// --------------------------------------------------
// MAIN TEST (constant field)
// --------------------------------------------------
int main()
{
    auto mesh = std::make_shared<TriangularMesh>("projects/Project-1/mesh_refined_2394.gri");
    mesh->writeGri("projects/Project-1/mesh_refined_2394");

    Eigen::MatrixXi I2E  = loadMatrix("projects/Project-1/mesh_refined_2394I2E.txt", 4).cast<int>();
    Eigen::MatrixXi B2E  = loadMatrix("projects/Project-1/mesh_refined_2394B2E.txt", 3).cast<int>();
    Eigen::MatrixXd In   = loadMatrix("projects/Project-1/mesh_refined_2394In.txt", 2);
    Eigen::MatrixXd Bn   = loadMatrix("projects/Project-1/mesh_refined_2394Bn.txt", 2);
    Eigen::VectorXd Area = loadVector("projects/Project-1/mesh_refined_2394Area.txt");

    std::cout << "I2E rows: " << I2E.rows() << "\n";
    std::cout << "In rows:  " << In.rows()  << "\n";
    std::cout << "B2E rows: " << B2E.rows() << "\n";
    std::cout << "Bn rows:  " << Bn.rows()  << "\n";

    // ---- bcNames (so bc is correct length) ----
    std::vector<std::string> bcNames;
    for (const auto& face : mesh->getFaces()) {
        if (!face.isBoundaryFace()) break;
        if (face.isPeriodicFace()) continue;
        if (std::find(bcNames.begin(), bcNames.end(), face._title) == bcNames.end())
            bcNames.push_back(face._title);
    }

    std::cout << "bcNames.size() = " << bcNames.size() << "\n";
    for (auto& s : bcNames) std::cout << "  " << s << "\n";

    // Use CopyBC (Ub = UL) for consistency in tests
    auto bc0 = std::make_shared<ConstantBC>(Eigen::Vector4d::Ones());
    std::vector<std::shared_ptr<BoundaryCondition>> bc(bcNames.size(), bc0);


    // ---- constant field ----
    StateMesh states(mesh, bc);
    states.state(0).fill(1.0);
    states.state(1).fill(1.0);
    states.state(2).fill(1.0);
    states.state(3).fill(1.0);

    const int nElem = states.cellCount();

    // ------------------ Compute FVGradient (WalkGrad) ------------------
    WalkGrad walkGrad;
    GradVec grads = walkGrad.computeGradient(I2E, B2E, In, Bn, Area, states);


    // // ---- FVGradient ----
    // HybridWalkPNGrad hybridWalkPNGrad;
    // auto grads = hybridWalkPNGrad.computeGrad_wpn(I2E, B2E, In, Bn, Area, states);

    StateMatrix gradX = StateMatrix::Zero(4, nElem);
    StateMatrix gradY = StateMatrix::Zero(4, nElem);
    for (int e = 0; e < nElem; ++e) {
        gradX.col(e) = grads[e].col(0);
        gradY.col(e) = grads[e].col(1);
    }

    StateMatrix Uplot = states.matrix();
    writeVTK(*mesh, Uplot, gradX, gradY, "gradient_output_constant.vtk");
    std::cout << "VTK file written: gradient_output_constant.vtk\n";

    // ---- verify ----
    double tol = 1e-12;
    bool pass = true;
    for (int e = 0; e < nElem; ++e)
        for (int v = 0; v < 4; ++v)
            if (std::abs(gradX(v,e)) > tol || std::abs(gradY(v,e)) > tol) {
                std::cout << "Nonzero grad at elem " << e << " var " << v
                          << " gradX=" << gradX(v,e) << " gradY=" << gradY(v,e) << "\n";
                pass = false;
            }

    if (pass) std::cout << "Gradient constant-field test PASSED\n";
    else      std::cout << "Gradient constant-field test FAILED\n";

    return pass ? 0 : 1;
}
