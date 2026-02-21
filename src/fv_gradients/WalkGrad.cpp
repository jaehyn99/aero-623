// Implementation of the gradient function using walk over 3 edges for both interior cells and boundary cells

# include "fv_gradients/WalkGrad.h"
# include "mesh/TriangularMesh.h"
# include "mesh/StateMesh.h"
# include "boundary_condition/BoundaryCondition.h"

// not sure if  i should keep this ======================
#include <Eigen/StdVector>
// ======================================================


# include <Eigen/Dense>
# include <Eigen/StdVector>

// ======================================================
// Dubug Includes
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <iostream>
// ======================================================

using State = Eigen::Vector4d;


// ======================================================
// GRADIENT CONSTRUCTION FUNCTION
// ======================================================
std::vector<Eigen::Matrix<double,4,2>> WalkGrad::computeGradient(const Eigen::MatrixXi& I2E,
                                                                const Eigen::MatrixXi& B2E,
                                                                const Eigen::MatrixXd& In,
                                                                const Eigen::MatrixXd& Bn,
                                                                const Eigen::VectorXd& Area,
                                                                const StateMesh& U) const
{

    // open once per call (overwrite each run)
    std::ofstream dbg("periodic_edges_debug.txt", std::ios::trunc);
    dbg << std::setprecision(17);
    dbg << "Periodic edge debug log\n";

    Eigen::MatrixXd gradX = Eigen::MatrixXd::Zero(U.stateCount(), U.cellCount());
    Eigen::MatrixXd gradY = Eigen::MatrixXd::Zero(U.stateCount(), U.cellCount());

    auto mesh = U.mesh();

    std::vector<bool> visited(mesh->numFaces(), false);

    std::vector<std::string> bcNames;
    // Loops through the boundary edges and add their names
    for (const auto& face: mesh->getFaces()){
        if (!face.isBoundaryFace()) break;
        if (face.isPeriodicFace()) continue;
        auto it = std::find(bcNames.cbegin(), bcNames.cend(), face._title);
        if (it == bcNames.cend()) bcNames.push_back(face._title);
    }

    // ===============================
    // Interior Faces
    // ===============================
    for (int f_i = 0; f_i<I2E.rows(); ++f_i)
    {
        // Get left and right cell indices for the edge
        int elemL = I2E(f_i,0)-1;
        int elemR = I2E(f_i,2)-1;
        int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
        int faceID = mesh->elem(elemL)._faceID[lfL];  // global face index

        const auto& face = mesh->face(faceID);

         // DEBUG BLOCK
        Eigen::Vector2d p0 = mesh->node(face._pointID[0]);
        Eigen::Vector2d p1 = mesh->node(face._pointID[1]);

        Eigen::Vector2d t = p1 - p0;
        double Lgeom = t.norm();

        // a unit geometric normal (one of the two possible signs)
        Eigen::Vector2d n_geom(t(1), -t(0));
        if (Lgeom > 0) n_geom /= Lgeom;

        // normalize the file normal for comparison
        Eigen::Vector2d n_file = In.row(f_i).transpose();
        double nf = n_file.norm();
        if (nf > 0) n_file /= nf;

        double absdot = std::abs(n_geom.dot(n_file));

        static int bad = 0;
        if (absdot < 0.95 && bad < 20) {
            std::cout << "[BAD MAP] f_i=" << f_i
                    << " faceID=" << faceID
                    << " abs(dot)=" << absdot
                    << "  n_file=" << n_file.transpose()
                    << "  n_geom=" << n_geom.transpose()
                    << "  Lgeom=" << Lgeom
                    << "  In_norm_raw=" << In.row(f_i).norm()
                    << "\n";
            bad++;
        }
        // END DEBUG BLOCK

        // DEBUG BLOCK
        static int p = 0;
        if (p < 10) {
            // compute n_geom for interior face exactly the same way
            std::cout << "f_i=" << f_i
                    << " abs(dot)=" << absdot
                    << " In=" << In.row(f_i)
                    << " n_geom=" << n_geom.transpose()
                    << "\n";
            p++;
        }
        // END DEBUG BLOCK
        
        
        if (visited[faceID]) {
            continue; // Skip if this face has already been visited (handles periodic pairs)
        }
        visited[faceID] = true;
        if (face.isPeriodicFace())
            visited[face._periodicFaceID] = true; // Mark the periodic counterpart as visited to avoid double counting

        // Get edge normal vector and length
        Eigen::Vector2d normal = In.row(f_i);
        double edge_length = mesh->length(faceID);

        if (face.isPeriodicFace()){
            if (face._elemID[0] > face._periodicElemID){
                normal *= -1;
            }
        }

        // Get cell center state values for left and right cells
        /** Instead of storing a copy of the state vector, we can directly reference the elements */
        State UL = U.cell(elemL);
        State UR = U.cell(elemR);

        // Compute state value at edge midpoint as average of left and right cell center values
        State Uhat = 0.5 * (UL + UR);

        // Compute gradient contribution from edge and add to left cell gradient and subtract from right cell gradient
        gradX.col(elemL) += Uhat * normal(0) * edge_length;
        gradY.col(elemL) += Uhat * normal(1) * edge_length;
        gradX.col(elemR) -= Uhat * normal(0) * edge_length; // opposite sign for right cell
        gradY.col(elemR) -= Uhat * normal(1) * edge_length; // opposite sign for right cell

        // Check if edge is periodic
        const bool isPeriodic = (face._periodicFaceID != -1);   // periodic interior connection

        // Dump if periodic
        if (isPeriodic) {
            const Eigen::Vector2d cL = mesh->centroid(elemL);
            const Eigen::Vector2d cR = mesh->centroid(elemR);

            dbg << "\n=== PERIODIC EDGE HIT (interior loop) ===\n";
            dbg << "edge(f_i) = " << f_i << " / " << I2E.rows() << "\n";
            dbg << "faceID    = " << faceID << "\n";

            dbg << "L: elemID = " << elemL << "   centroid = " << cL.transpose() << "\n";
            dbg << "R: elemID = " << elemR << "   centroid = " << cR.transpose() << "\n";

            dbg << "lfL       = " << lfL << "\n";
            dbg << "I2E row    = " << I2E(f_i,0) << " " << I2E(f_i,1) << " "
                            << I2E(f_i,2) << " " << I2E(f_i,3) << "\n";

            dbg << "normal    = " << normal.transpose() << "\n";
            dbg << "length    = " << edge_length << "\n";

            dbg << "face pts  = " << face._pointID[0] << " " << face._pointID[1] << "\n";
            dbg << "elemID    = " << face._elemID[0] << " " << face._elemID[1] << "\n";
            dbg << "periodicFaceID = " << face._periodicFaceID << "\n";
            dbg << "periodicElemID = " << face._periodicElemID << "\n";
            dbg << "title     = " << face._title << "\n";
            dbg << "nf        = " << face._nf << "\n";

            dbg << "UL        = " << UL.transpose() << "\n";
            dbg << "UR        = " << UR.transpose() << "\n";

            dbg << "gradX(L) before = " << gradX.col(elemL).transpose() << "\n";
            dbg << "gradY(L) before = " << gradY.col(elemL).transpose() << "\n";
            dbg << "gradX(R) before = " << gradX.col(elemR).transpose() << "\n";
            dbg << "gradY(R) before = " << gradY.col(elemR).transpose() << "\n";
            dbg << "========================================\n";
            dbg.flush(); // optional
        }
    }

    // ===============================
    // Boundary Faces
    // ===============================
    static int hit = 0;
    for (int f_b = 0; f_b < B2E.rows(); ++f_b)
    {
        // DEBUG BLOCK
        
        if (hit == 0) std::cout << "Entered boundary loop\n";
        hit++;
        // END DEBUG BLOCK

        // Get left cell index for the edge (there is no right cell since this is a boundary edge)
        int elem   = B2E(f_b,0)-1;
        int lf     = B2E(f_b,1) - 1;
        int faceID = mesh->elem(elem)._faceID[lf];
        int bgroup = B2E(f_b,2);

        // DEBUG BLOCK
        static int printed = 0;
        if (printed < 300) {
            const auto& face = mesh->face(faceID);
            auto it = std::find(bcNames.cbegin(), bcNames.cend(), face._title);
            std::size_t boundaryID = std::size_t(it - bcNames.cbegin());

            std::cout << "f_b=" << f_b
                    << " title=" << face._title
                    << " bgroup=" << bgroup
                    << " boundaryID=" << boundaryID
                    << "\n";
            printed++;
        }
        // END DEBUG BLOCK

       


        // Get edge normal vector and length
        Eigen::Vector2d normal = Bn.row(f_b);
        double edge_length = mesh->length(faceID);

        // Get cell center state values for left cell
        State UL = U.cell(elem);

        // Use boundary state based on bgroup (CHANGE BGROUP NUMBERS TO MATCH YOUR MESH)
        const auto& face = mesh->face(faceID);

        // --- DEBUG: ensure Bn points outward for this boundary face
        {

            // Face endpoints and midpoint
            Eigen::Vector2d p0 = mesh->node(face._pointID[0]);
            Eigen::Vector2d p1 = mesh->node(face._pointID[1]);
            Eigen::Vector2d xf = 0.5 * (p0 + p1);

            // Cell centroid
            Eigen::Vector2d xc = mesh->centroid(elem);

            // A geometric unit normal (one of the two signs)
            Eigen::Vector2d t = p1 - p0;
            double L = t.norm();
            Eigen::Vector2d n_geom(t(1), -t(0));
            if (L > 0.0) n_geom /= L;

            // Choose the outward sign: outward points from centroid toward face
            if ((xf - xc).dot(n_geom) < 0.0) n_geom *= -1.0;

            // File normal (unit, per your stats)
            Eigen::Vector2d n_file = Bn.row(f_b).transpose();

            // If dot < 0, Bn is pointing inward relative to this cell
            double dot = n_geom.dot(n_file);

            static int printed = 0;
            if (printed < 10) {
                std::cout << "f_b=" << f_b
                        << " bgroup=" << bgroup
                        << " title=" << face._title
                        << " dot=" << dot
                        << " n_file=" << n_file.transpose()
                        << " n_out=" << n_geom.transpose()
                        << "\n";
                printed++;
            }

            // OPTIONAL FIX (try this): force outward normal in your computation
            // Eigen::Vector2d normal = n_file;
            // if (dot < 0.0) normal *= -1.0;
        }
        // --- END DEBUG

        std::size_t boundaryID = std::find(bcNames.cbegin(), bcNames.cend(), face._title) - bcNames.cbegin();
        State Ub = U.bc(boundaryID)->computeBoundaryState(UL, normal);
        // State Ub;
        // if (bgroup == 3)
        //     Ub = inletState_(UL, gamma_, normal);
        // else if (bgroup == 7)
        //     Ub = outletState_(UL, gamma_, normal);
        // else if (bgroup == 1 || bgroup == 5)
        //     Ub = wallState_(UL, gamma_, normal);
        // else
        //     throw std::runtime_error("Unknown boundary group in computeGradient");

        // Compute gradient contribution from edge and add to left cell gradient
        gradX.col(elem) += Ub * normal(0) * edge_length;
        gradY.col(elem) += Ub * normal(1) * edge_length;
    }

    // DEBUG BLOCK
    std::cout << "Boundary loop iterations = " << hit << "\n";
    // END DEBUG BLOCK

    // ===============================
    // Divide gradient by Area
    // ===============================
    for (int e = 0; e < U.cellCount(); ++e)
    {
        gradX.col(e) /= Area(e);
        gradY.col(e) /= Area(e);
    }

    std::vector<Eigen::Matrix<double,4,2>> Lgrad;
    Lgrad.resize(U.cellCount());

    for (int e = 0; e < U.cellCount(); ++e)
    {
        Lgrad[e].col(0) = gradX.col(e);   // dx
        Lgrad[e].col(1) = gradY.col(e);   // dy
    }

    return Lgrad;
}



// # include "solver/SecondorderEuler.h"
// # include "mesh/TriangularMesh.h"

// # include <Eigen/Dense>

// using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
// using State = Eigen::Vector4d;

// // ======================================================
// // Dubug Includes
// #include <fstream>
// #include <iomanip>
// #include <cstdlib>
// // ======================================================

// // Outer Main Function Call
// StateMatrix SecondOrderEuler::computeResidual(
//     TriangularMesh& mesh,
//     const Eigen::MatrixXi& I2E,
//     const Eigen::MatrixXi& B2E,
//     const Eigen::MatrixXd& In,
//     const Eigen::MatrixXd& Bn,
//     const Eigen::VectorXd& Area,
//     const StateMatrix& U) const
// {
//     const int nElem = U.cols();
//     const int nVar  = U.rows();

//     StateMatrix gradX = StateMatrix::Zero(nVar, nElem);
//     StateMatrix gradY = StateMatrix::Zero(nVar, nElem);

//     computeGradient(mesh, I2E, B2E, In, Bn, Area, U, gradX, gradY);

//     return computeResidualFromGradient(mesh, I2E, B2E, In, Bn, Area,
//                                        U, gradX, gradY);
// }


// // ======================================================
// // GRADIENT CONSTRUCTION FUNCTION
// // ======================================================
// void SecondOrderEuler::computeGradient(
//     TriangularMesh& mesh,
//     const Eigen::MatrixXi& I2E,
//     const Eigen::MatrixXi& B2E,
//     const Eigen::MatrixXd& In,
//     const Eigen::MatrixXd& Bn,
//     const Eigen::VectorXd& Area,
//     const StateMatrix& U,
//     StateMatrix& gradX,
//     StateMatrix& gradY) const
// {
//     // open once per call (overwrite each run)
//     std::ofstream dbg("periodic_edges_debug.txt", std::ios::trunc);
//     dbg << std::setprecision(17);
//     dbg << "Periodic edge debug log\n";

//     std::vector<bool> visited(mesh.numFaces(), false);
//     // ===============================
//     // ===============================
//     // Interior Faces
//     // ===============================
//     for (int f_i = 0; f_i<I2E.rows(); ++f_i)
//     {
//         // Get left and right cell indices for the edge
//         int elemL = I2E(f_i,0)-1;
//         int elemR = I2E(f_i,2)-1;
//         int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
//         int faceID = mesh.elem(elemL)._faceID[lfL];  // global face index
        
//         const auto& face = mesh.face(faceID);
//         if (visited[faceID]) {
//             continue; // Skip if this face has already been visited (handles periodic pairs)
//         }
//         visited[faceID] = true;
//         visited[face._periodicFaceID] = true; // Mark the periodic counterpart as visited to avoid double counting

//         // Get edge normal vector and length
//         Eigen::Vector2d normal = In.row(f_i);
//         double edge_length = mesh.length(faceID);

//         if (face.isPeriodicFace()){
//             if (face._elemID[0] > face._periodicElemID){
//                 normal *= -1;
//             }
//         }

//         // Get cell center state values for left and right cells
//         /** Instead of storing a copy of the state vector, we can directly reference the elements */
//         const State& UL = U.col(elemL);
//         const State& UR = U.col(elemR);

//         // Compute state value at edge midpoint as average of left and right cell center values
//         State Uhat = 0.5 * (UL + UR);

//         // Compute gradient contribution from edge and add to left cell gradient and subtract from right cell gradient
//         gradX.col(elemL) += Uhat * normal(0) * edge_length;
//         gradY.col(elemL) += Uhat * normal(1) * edge_length;
//         gradX.col(elemR) -= Uhat * normal(0) * edge_length; // opposite sign for right cell
//         gradY.col(elemR) -= Uhat * normal(1) * edge_length; // opposite sign for right cell

//         // Check if edge is periodic
//         const bool isPeriodic = (face._periodicFaceID != -1);   // periodic interior connection

//         // Dump if periodic
//         if (isPeriodic) {
//             const Eigen::Vector2d cL = mesh.centroid(elemL);
//             const Eigen::Vector2d cR = mesh.centroid(elemR);

//             dbg << "\n=== PERIODIC EDGE HIT (interior loop) ===\n";
//             dbg << "edge(f_i) = " << f_i << " / " << I2E.rows() << "\n";
//             dbg << "faceID    = " << faceID << "\n";

//             dbg << "L: elemID = " << elemL << "   centroid = " << cL.transpose() << "\n";
//             dbg << "R: elemID = " << elemR << "   centroid = " << cR.transpose() << "\n";

//             dbg << "lfL       = " << lfL << "\n";
//             dbg << "I2E row    = " << I2E(f_i,0) << " " << I2E(f_i,1) << " "
//                             << I2E(f_i,2) << " " << I2E(f_i,3) << "\n";

//             dbg << "normal    = " << normal.transpose() << "\n";
//             dbg << "length    = " << edge_length << "\n";

//             dbg << "face pts  = " << face._pointID[0] << " " << face._pointID[1] << "\n";
//             dbg << "elemID    = " << face._elemID[0] << " " << face._elemID[1] << "\n";
//             dbg << "periodicFaceID = " << face._periodicFaceID << "\n";
//             dbg << "periodicElemID = " << face._periodicElemID << "\n";
//             dbg << "title     = " << face._title << "\n";
//             dbg << "nf        = " << face._nf << "\n";

//             dbg << "UL        = " << UL.transpose() << "\n";
//             dbg << "UR        = " << UR.transpose() << "\n";

//             dbg << "gradX(L) before = " << gradX.col(elemL).transpose() << "\n";
//             dbg << "gradY(L) before = " << gradY.col(elemL).transpose() << "\n";
//             dbg << "gradX(R) before = " << gradX.col(elemR).transpose() << "\n";
//             dbg << "gradY(R) before = " << gradY.col(elemR).transpose() << "\n";
//             dbg << "========================================\n";
//             dbg.flush(); // optional
//         }
//     }

//     // ===============================
//     // Boundary Faces
//     // ===============================
//     for (int f_b = 0; f_b < B2E.rows(); ++f_b)
//     {
//         // Get left cell index for the edge (there is no right cell since this is a boundary edge)
//         int elem   = B2E(f_b,0)-1;
//         int lf     = B2E(f_b,1) - 1;
//         int faceID = mesh.elem(elem)._faceID[lf];
//         int bgroup = B2E(f_b,2);

//         // Get edge normal vector and length
//         Eigen::Vector2d normal = Bn.row(f_b);
//         double edge_length = mesh.length(faceID);

//         // Get cell center state values for left cell
//         const State& UL = U.col(elem);

//         // Use boundary state based on bgroup (CHANGE BGROUP NUMBERS TO MATCH YOUR MESH)
//         State Ub;
//         if (bgroup == 3)
//             Ub = inletState_(UL, gamma_, normal);
//         else if (bgroup == 7)
//             Ub = outletState_(UL, gamma_, normal);
//         else if (bgroup == 1 || bgroup == 5)
//             Ub = wallState_(UL, gamma_, normal);
//         else
//             throw std::runtime_error("Unknown boundary group in computeGradient");

//         // Compute gradient contribution from edge and add to left cell gradient
//         gradX.col(elem) += Ub * normal(0) * edge_length;
//         gradY.col(elem) += Ub * normal(1) * edge_length;
//     }

//     // ===============================
//     // Divide gradient by Area
//     // ===============================
//     for (int e = 0; e < U.cols(); ++e)
//     {
//         gradX.col(e) /= Area(e);
//         gradY.col(e) /= Area(e);
//     }
// }
// ======================================================
// END OF GRADIENT CONSTRUCTION FUNCTION
// ======================================================


// // ======================================================
// // RESIDUAL CONSTRUCTION FUNCTION
// // ======================================================
// StateMatrix SecondOrderEuler::computeResidualFromGradient(
//     TriangularMesh& mesh,
//     const Eigen::MatrixXi& I2E,
//     const Eigen::MatrixXi& B2E,
//     const Eigen::MatrixXd& In,
//     const Eigen::MatrixXd& Bn,
//     const Eigen::VectorXd& Area,
//     const StateMatrix& U,
//     const StateMatrix& gradX,
//     const StateMatrix& gradY) const
// {
//     std::vector<bool> visited(mesh.numFaces(), false);

//     StateMatrix R = StateMatrix::Zero(U.rows(), U.cols());
//     // ===============================
//     // Interior Faces
//     // ===============================
//     for (int f_i = 0; f_i < I2E.rows(); ++f_i)
//     {
//         // Get left and right cell indices for the edge
//         int elemL = I2E(f_i,0)-1;
//         int elemR = I2E(f_i,2)-1;
//         int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
//         int faceID = mesh.elem(elemL)._faceID[lfL];  // global face index

//         const auto& face = mesh.face(faceID);
//         if (visited[faceID]) {
//             continue; // Skip if this face has already been visited (handles periodic pairs)
//         }
//         visited[faceID] = true;
//         visited[face._periodicFaceID] = true; // Mark the periodic counterpart as visited to avoid double counting

//         // Get edge normal vector, length and midpoint coordinates
//         Eigen::Vector2d normal = In.row(f_i);
//         double edge_length = mesh.length(faceID);
//         // const auto& face = mesh.face(faceID);
//         // Eigen::Vector2d midpoint = 0.5 * (mesh.node(face._pointID[0]) + mesh.node(face._pointID[1]));

//         // Get state vector at the edge midpoint for left face
//         int faceIDL = mesh.elem(elemL)._faceID[lfL];
//         const auto& faceL = mesh.face(faceIDL);
//         Eigen::Vector2d midpointL = 0.5 * (mesh.node(faceL._pointID[0]) + mesh.node(faceL._pointID[1]));
//         const State& UL = U.col(elemL);
//         Eigen::Vector2d centroidL = mesh.centroid(elemL);
//         Eigen::Vector2d dxL = midpointL - centroidL;
//         State UhatL = UL + gradX.col(elemL) * dxL(0) + gradY.col(elemL) * dxL(1);

//         // Get state vector at the edge midpoint for right face
//         int lfR = I2E(f_i,3) - 1; // local face on elemR
//         int faceIDR = mesh.elem(elemR)._faceID[lfR];
//         const auto& faceR = mesh.face(faceIDR);
//         Eigen::Vector2d midpointR = 0.5 * (mesh.node(faceR._pointID[0]) + mesh.node(faceR._pointID[1]));
//         const State& UR = U.col(elemR);
//         Eigen::Vector2d centroidR = mesh.centroid(elemR);
//         Eigen::Vector2d dxR = midpointR - centroidR;
//         State UhatR = UR + gradX.col(elemR) * dxR(0) + gradY.col(elemR) * dxR(1);

//         // use flux function to get the flux at the edge midpoint using uL and uR
//         Eigen::Vector4d flux = numFlux_(UhatL, UhatR, gamma_, normal);
//         flux *= edge_length;

//         // add the flux contribution to the residual of the left cell and subtract it from the residual of the right cell. 
//             // When this is done over all the edges this will construct the full residual over all the cells in the domain
//         R.col(elemL) += flux;
//         R.col(elemR) -= flux;   // opposite sign

//         // compute the wave speed at the edge and update the maximum wave speed for the left and right cells if necessary
        
//     }

//     // ===============================
//     // Boundary Faces
//     // ===============================
//     for (int f_b = 0; f_b < B2E.rows(); ++f_b)
//     {
//         // Get left cell index for the edge (there is no right cell since this is a boundary edge)
//         int elem   = B2E(f_b,0)-1;
//         int lf     = B2E(f_b,1) - 1;
//         int faceID = mesh.elem(elem)._faceID[lf];
//         int bgroup = B2E(f_b,2);

//         // Get edge normal vector and length
//         Eigen::Vector2d normal = Bn.row(f_b);
//         double edge_length = mesh.length(faceID);

//         // Get cell center state values for left cell
//         const State& UL = U.col(elem);

//         // Based on boundary conditions, choose flux (CHANGE BGROUP NUMBERS TO MATCH YOUR MESH)
//         if (bgroup == 3) { // inlet
//             // compute the flux at the boundary using the inlet flux function and the interior state at the cell center and the boundary state
//             Eigen::Vector4d flux = inletFlux_(UL, gamma_, normal);
//             flux *= edge_length;
//             R.col(elem) += flux; // add contribution to residual of left cell
//         }
//         else if (bgroup == 7) { // outlet
//             Eigen::Vector4d flux = outletFlux_(UL, gamma_, normal);
//             flux *= edge_length;
//             R.col(elem) += flux; // add contribution to residual of left cell
//         }
//         else if (bgroup == 1 || bgroup == 5) { // wall
//             Eigen::Vector4d flux = wallFlux_(UL, gamma_, normal);
//             flux *= edge_length;
//             R.col(elem) += flux; // add contribution to residual of left cell
//         }
//         else
//             throw std::runtime_error("Unknown boundary group in computeResidualFromGradient");

//     }

//     return R;
// }

// } // namespace solver