#include <iostream>
#include <cassert>
#include <fstream>
#include <Eigen/Dense>

#include "solver/SecondorderEuler.h"
#include "mesh/TriangularMesh.h"
#include "CopyInteriorState.h"

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
              const Eigen::MatrixXd& U,
              const Eigen::MatrixXd& gradX,
              const Eigen::MatrixXd& gradY,
              const std::string& filename)
{
    std::ofstream out(filename);

    int nNodes = mesh.numNodes();
    int nElems = mesh.numElems();

    out << "# vtk DataFile Version 3.0\n";
    out << "State and Gradient Output\n";
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

    // ---- Write State Variables ----
    for (int v = 0; v < 4; ++v)
    {
        out << "SCALARS U_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e)
            out << U(v,e) << "\n";
    }

    // ---- Write gradX ----
    for (int v = 0; v < 4; ++v)
    {
        out << "SCALARS gradX_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e)
            out << gradX(v,e) << "\n";
    }

    // ---- Write gradY ----
    for (int v = 0; v < 4; ++v)
    {
        out << "SCALARS gradY_var" << v << " double\n";
        out << "LOOKUP_TABLE default\n";
        for (int e = 0; e < nElems; ++e)
            out << gradY(v,e) << "\n";
    }

    out.close();
}


// --------------------------------------------------
// Dummy fluxes
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
    TriangularMesh mesh("projects/Project-1/meshGlobalRefined1.gri");
    mesh.writeGri("projects/Project-1/meshGlobalRefined1_output.gri"); // Write out to verify correct loading

    // ------------------ Load Connectivity ------------------
    Eigen::MatrixXi I2E = loadMatrix("projects/Project-1/meshGlobalRefined1_outputI2E.txt", 4).cast<int>();
    Eigen::MatrixXi B2E = loadMatrix("projects/Project-1/meshGlobalRefined1_outputB2E.txt", 3).cast<int>();
    Eigen::MatrixXd In  = - loadMatrix("projects/Project-1/meshGlobalRefined1_outputIn.txt", 2);
    Eigen::MatrixXd Bn  = - loadMatrix("projects/Project-1/meshGlobalRefined1_outputBn.txt", 2);
    Eigen::VectorXd Area = loadVector("projects/Project-1/meshGlobalRefined1_outputArea.txt");

    // ------------------ Solver Setup ------------------
    DummyFlux numFlux;
    DummyBoundaryFlux inletFlux;
    DummyBoundaryFlux outletFlux;
    DummyBoundaryFlux wallFlux;

    CopyInteriorState inletState;
    CopyInteriorState outletState;
    CopyInteriorState wallState;

    solver::SecondOrderEuler solver(numFlux, inletFlux, outletFlux, wallFlux, inletState, outletState, wallState, 1.4);

    int nElem = Area.size();

    // --------------------------------------------------
    // Linear X-Variation Test
    // U(x,y) = slope * x
    // --------------------------------------------------

    StateMatrix U(4, nElem);
    StateMatrix gradX = StateMatrix::Zero(4, nElem);
    StateMatrix gradY = StateMatrix::Zero(4, nElem);

    Eigen::Vector4d slope;
    slope << 1.0, 2.0, -0.5, 3.0;

    for (int e = 0; e < nElem; ++e)
    {
        Eigen::Vector2d centroid = mesh.centroid(e);
        double x = centroid(0);

        for (int v = 0; v < 4; ++v)
            U(v,e) = slope(v) * x;
    }

    // // Walk over 3 edges test
    // solver.test_ComputeGradient(mesh, I2E, B2E, In, Bn, Area, U, gradX, gradY);

    // Plane normals test
    solver.test_ComputeGradient_pn(mesh, I2E, B2E, In, Bn, Area, U, gradX, gradY);

    // ------------------ Write VTK ------------------
    writeVTK(mesh, U, gradX, gradY, "gradient_output_linear.vtk");

    std::cout << "VTK file written: gradient_output_linear.vtk\n";

    // --------------------------------------------------
    // Verification
    // --------------------------------------------------

    double tol = 1e-10;
    bool pass = true;

    double sumSqX = 0.0;
    double sumSqY = 0.0;
    int count = 0;

    for (int e = 0; e < nElem; ++e)
    {
        for (int v = 0; v < 4; ++v)
        {
            double errX = gradX(v,e) - slope(v);
            double errY = gradY(v,e);   // expected 0

            sumSqX += errX * errX;
            sumSqY += errY * errY;
            count++;

            if (std::abs(errX) > tol)
            {
                std::cout << "gradX mismatch at elem "
                        << e << " var " << v
                        << " Expected: " << slope(v)
                        << " Got: " << gradX(v,e) << "\n";
                pass = false;
            }

            if (std::abs(errY) > tol)
            {
                std::cout << "gradY nonzero at elem "
                        << e << " var " << v
                        << " Value: " << gradY(v,e) << "\n";
                pass = false;
            }
        }
    }

    double rmsX = std::sqrt(sumSqX / count);
    double rmsY = std::sqrt(sumSqY / count);

    std::cout << "\nRMS gradX error: " << rmsX << "\n";
    std::cout << "RMS gradY error: " << rmsY << "\n";



    if (pass)
        std::cout << "Linear X-gradient test PASSED\n";
    else
        std::cout << "Linear X-gradient test FAILED\n";

    return 0;
}
