#include <iostream>
#include <cassert>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <cmath>
#include <vector>

#include "mesh/TriangularMesh.h"
#include "CopyBC.h"
#include "mesh/StateMesh.h"
#include "fv_gradients/WalkGrad.h"
#include "fv_gradients/hybridWalkPNGrad.h"
#include "fv_gradients/Gradients.h"
#include "boundary_condition/InletOutletBC.h"
#include "boundary_condition/InviscidWallBC.h"
#include "boundary_condition/OutletBC.h"

using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using GradCell = Eigen::Matrix<double,4,2>;
using GradVec  = std::vector<Eigen::Matrix<double,4,2>>;

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
              const StateMatrix& U,
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
// MAIN TEST
// --------------------------------------------------
int main()
{

    // ------------------ Load Mesh ------------------
    auto mesh = std::make_shared<TriangularMesh>("projects/Project-1/mesh_refined_2394.gri");
    mesh->writeGri("projects/Project-1/mesh_refined_2394");

    // ------------------ Load Connectivity ------------------
    Eigen::MatrixXi I2E  = loadMatrix("projects/Project-1/mesh_refined_2394I2E.txt", 4).cast<int>();
    Eigen::MatrixXi B2E  = loadMatrix("projects/Project-1/mesh_refined_2394B2E.txt", 3).cast<int>();
    Eigen::MatrixXd In   = loadMatrix("projects/Project-1/mesh_refined_2394In.txt", 2);
    Eigen::MatrixXd Bn   = loadMatrix("projects/Project-1/mesh_refined_2394Bn.txt", 2);
    Eigen::VectorXd Area = loadVector("projects/Project-1/mesh_refined_2394Area.txt");

    // DEBUG BLOCK
    auto stats = [](const Eigen::MatrixXd& N, const std::string& name){
        double mn = 1e300, mx = -1e300, mean = 0.0;
        for (int i=0;i<N.rows();++i){
            double s = N.row(i).norm();
            mn = std::min(mn,s); mx = std::max(mx,s); mean += s;
        }
        mean /= N.rows();
        std::cout << name << " row-norm stats: min=" << mn
                << " mean=" << mean
                << " max=" << mx << "\n";
    };
    // END DEBUG BLOCK

    stats(In,"In");
    stats(Bn,"Bn");


    std::cout << "I2E rows: " << I2E.rows() << "\n";
    std::cout << "In rows:  " << In.rows()  << "\n";
    std::cout << "B2E rows: " << B2E.rows() << "\n";
    std::cout << "Bn rows:  " << Bn.rows()  << "\n";

    // Build bcs
    auto copy = std::make_shared<CopyBC>();



    // Debug ==============================
    std::vector<std::string> bcNames;
    for (const auto& face : mesh->getFaces()) {
        if (!face.isBoundaryFace()) break;
        if (face.isPeriodicFace()) continue;
        if (std::find(bcNames.begin(), bcNames.end(), face._title) == bcNames.end())
            bcNames.push_back(face._title);
    }

    std::cout << "bcNames.size() = " << bcNames.size() << "\n";
    for (auto& s : bcNames) std::cout << "  " << s << "\n";
    // =====================================

    std::vector<std::shared_ptr<BoundaryCondition>> bc(bcNames.size(), copy);

    // // Inlet conditions
    // double gamma = 1.4;
    // double rho0 = 1;
    // double a0 = 1;
    // double p0 = rho0*a0*a0/gamma;
    // double alpha = 50*3.14/180;
    // double pout = 0.7*p0;
    // double M = 0.1;

    // std::shared_ptr<InletOutletBC> inlet = std::make_shared<InletOutletBC>(rho0, a0, alpha, pout, gamma);
    // std::shared_ptr<BoundaryCondition> wall = std::make_shared<InviscidWallBC>(gamma);
    // std::shared_ptr<BoundaryCondition> outlet = std::make_shared<OutletBC>(pout, gamma);
    // std::vector<std::shared_ptr<BoundaryCondition>> bc{wall, inlet, wall, outlet};
    // Eigen::Vector4d u{1, 1, 1, 5};
    // Eigen::Vector2d n{1, 0};
    // std::cout << bc[0]->computeFlux(u, n) << std::endl;

    // std::vector<std::shared_ptr<BoundaryCondition>> bc{copy, copy, copy, copy};

    // ------------------ Linear State Field Test ------------------
    const int nElem = int(Area.size());

    // Initialize the state mesh (same as before)
    StateMesh states(mesh, bc);

    // choose linear field coefficients for each state k
    Eigen::Vector4d ax; ax << 1.0, 2.0, -0.5, 3.0;      // dU/dx

    for (int e = 0; e < states.cellCount(); ++e)
    {
        Eigen::Vector2d xc = mesh->centroid(e);   // (x,y) of cell centroid
        double x = xc(0);                      // use x-coordinate for linear variation

        for (int k = 0; k < states.stateCount(); ++k)
            states(k,e) = ax(k)*x;
    }

    // ------------------ Solver Setup ------------------
    
    // ------------------ Compute gradients (WalkGrad) ------------------
    std::shared_ptr<Gradients> walkGrad = std::make_shared<WalkGrad>();
    states.setGradientMethod(walkGrad);
    GradVec grads = states.getGradients();

    // // ------------------ Compute gradients (Hybrid Walk Pn Grad) ------------------
    // std::shared_ptr<Gradients> hybridWalkPNGrad = std::make_shared<HybridWalkPNGrad>();
    // states.setGradientMethod(hybridWalkPNGrad);
    // GradVec grads = states.getGradients();

    // Convert to old matrix form for VTK + assertions
    StateMatrix gradX = StateMatrix::Zero(4, nElem);
    StateMatrix gradY = StateMatrix::Zero(4, nElem);
    for (int e = 0; e < nElem; ++e) {
        gradX.col(e) = grads[e].col(0);
        gradY.col(e) = grads[e].col(1);
    }
    StateMatrix Uplot = states.matrix();

    // ------------------ Write VTK ------------------
    writeVTK(*mesh, Uplot, gradX, gradY, "gradient_output_linear.vtk");

    std::cout << "VTK file written: gradient_output_linear.vtk\n";

    // --------------------------------------------------
    // Verification
    // --------------------------------------------------

    double tol = 1e-10;
    bool pass = true;

    double sumSqX = 0.0, sumSqY = 0.0;
    int count = 0;

    for (int e = 0; e < nElem; ++e)
    {
        for (int v = 0; v < 4; ++v)
        {
            double errX = gradX(v,e) - ax(v);

            sumSqX += errX * errX;
            count++;

            if (std::abs(errX) > tol)
            {
                std::cout << "gradX mismatch at elem "
                        << e << " var " << v
                        << " Expected: " << ax(v)
                        << " Got: " << gradX(v,e) << "\n";
                pass = false;
            }
        }
    }

    double rmsX = std::sqrt(sumSqX / count);

    std::cout << "\nRMS gradX error: " << rmsX << "\n";

    if (pass) std::cout << "Linear gradient test PASSED\n";
    else      std::cout << "Linear gradient test FAILED\n";


    return 0;
}
