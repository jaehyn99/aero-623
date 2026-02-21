# include "solver/SecondorderEuler.h"
# include "mesh/TriangularMesh.h"

# include <Eigen/Dense>

using GradMatrix = std::vector<Eigen::Matrix<double,4,2>>;
using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using State = Eigen::Vector4d;


// ======================================================
// RESIDUAL CONSTRUCTION FUNCTION
// ======================================================
Eigen::MatrixXd SecondOrderEuler::computeResidual(
    const Eigen::MatrixXi& I2E,
    const Eigen::MatrixXi& B2E,
    const Eigen::MatrixXd& In,
    const Eigen::MatrixXd& Bn,
    const Eigen::VectorXd& Area,
    const StateMesh& U,
    const GradMatrix& gradients) const
{
    std::vector<bool> visited(mesh.numFaces(), false);

    StateMatrix R = StateMatrix::Zero(U.rows(), U.cols());
    // ===============================
    // Interior Faces
    // ===============================
    for (int f_i = 0; f_i < I2E.rows(); ++f_i)
    {
        // Get left and right cell indices for the edge
        int elemL = I2E(f_i,0)-1;
        int elemR = I2E(f_i,2)-1;
        int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
        int faceID = mesh.elem(elemL)._faceID[lfL];  // global face index

        const auto& face = mesh.face(faceID);
        if (visited[faceID]) {
            continue; // Skip if this face has already been visited (handles periodic pairs)
        }
        visited[faceID] = true;
        visited[face._periodicFaceID] = true; // Mark the periodic counterpart as visited to avoid double counting

        // Get edge normal vector, length and midpoint coordinates
        Eigen::Vector2d normal = In.row(f_i);
        double edge_length = mesh.length(faceID);
        // const auto& face = mesh.face(faceID);
        // Eigen::Vector2d midpoint = 0.5 * (mesh.node(face._pointID[0]) + mesh.node(face._pointID[1]));

        // Get state vector at the edge midpoint for left face
        int faceIDL = mesh.elem(elemL)._faceID[lfL];
        const auto& faceL = mesh.face(faceIDL);
        Eigen::Vector2d midpointL = 0.5 * (mesh.node(faceL._pointID[0]) + mesh.node(faceL._pointID[1]));
        const State& UL = U.col(elemL);
        Eigen::Vector2d centroidL = mesh.centroid(elemL);
        Eigen::Vector2d dxL = midpointL - centroidL;
        State UhatL = UL + gradX.col(elemL) * dxL(0) + gradY.col(elemL) * dxL(1);

        // Get state vector at the edge midpoint for right face
        int lfR = I2E(f_i,3) - 1; // local face on elemR
        int faceIDR = mesh.elem(elemR)._faceID[lfR];
        const auto& faceR = mesh.face(faceIDR);
        Eigen::Vector2d midpointR = 0.5 * (mesh.node(faceR._pointID[0]) + mesh.node(faceR._pointID[1]));
        const State& UR = U.col(elemR);
        Eigen::Vector2d centroidR = mesh.centroid(elemR);
        Eigen::Vector2d dxR = midpointR - centroidR;
        State UhatR = UR + gradX.col(elemR) * dxR(0) + gradY.col(elemR) * dxR(1);

        // use flux function to get the flux at the edge midpoint using uL and uR
        Eigen::Vector4d flux = numFlux_(UhatL, UhatR, gamma_, normal);
        flux *= edge_length;

        // add the flux contribution to the residual of the left cell and subtract it from the residual of the right cell. 
            // When this is done over all the edges this will construct the full residual over all the cells in the domain
        R.col(elemL) += flux;
        R.col(elemR) -= flux;   // opposite sign

        // compute the wave speed at the edge and update the maximum wave speed for the left and right cells if necessary
        
    }

    // ===============================
    // Boundary Faces
    // ===============================
    for (int f_b = 0; f_b < B2E.rows(); ++f_b)
    {
        // Get left cell index for the edge (there is no right cell since this is a boundary edge)
        int elem   = B2E(f_b,0)-1;
        int lf     = B2E(f_b,1) - 1;
        int faceID = mesh.elem(elem)._faceID[lf];
        int bgroup = B2E(f_b,2);

        // Get edge normal vector and length
        Eigen::Vector2d normal = Bn.row(f_b);
        double edge_length = mesh.length(faceID);

        // Get cell center state values for left cell
        const State& UL = U.col(elem);

        // Based on boundary conditions, choose flux (CHANGE BGROUP NUMBERS TO MATCH YOUR MESH)
        if (bgroup == 3) { // inlet
            // compute the flux at the boundary using the inlet flux function and the interior state at the cell center and the boundary state
            Eigen::Vector4d flux = inletFlux_(UL, gamma_, normal);
            flux *= edge_length;
            R.col(elem) += flux; // add contribution to residual of left cell
        }
        else if (bgroup == 7) { // outlet
            Eigen::Vector4d flux = outletFlux_(UL, gamma_, normal);
            flux *= edge_length;
            R.col(elem) += flux; // add contribution to residual of left cell
        }
        else if (bgroup == 1 || bgroup == 5) { // wall
            Eigen::Vector4d flux = wallFlux_(UL, gamma_, normal);
            flux *= edge_length;
            R.col(elem) += flux; // add contribution to residual of left cell
        }
        else
            throw std::runtime_error("Unknown boundary group in computeResidualFromGradient");

    }

    return R;
}