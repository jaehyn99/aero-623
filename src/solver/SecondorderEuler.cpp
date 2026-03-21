// // WHAT I HAVE IMPLEMENTED HERE IS USING BOUNDARY STATE FOR GRADIENT CONSTRUCTION AT BOUNDARY CELLS

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


// namespace solver {

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
// // ======================================================
// // END OF GRADIENT CONSTRUCTION FUNCTION
// // ======================================================


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



// // StateMatrix
// // SecondOrderEuler::computeResidual(const TriangularMesh& mesh,
// //                                   const Eigen::MatrixXi& I2E,
// //                                   const Eigen::MatrixXi& B2E,
// //                                   const Eigen::MatrixXd& In,
// //                                   const Eigen::MatrixXd& Bn,
// //                                   const Eigen::VectorXd& Area,
// //                                   const StateMatrix& U) const
// // {

// //     // Can be replaced by other meshes
// //     ::TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");

// //     // Get number of elements and number of variables from the input state vector
// //     const int nElem = U.cols();
// //     const int nVar  = U.rows();

// //     // initialize the residual on all cells to zero
// //     StateMatrix R = StateMatrix::Zero(nVar, nElem);

// //     // initialize the gradient on all cells to zero
// //     StateMatrix gradX = StateMatrix::Zero(nVar, nElem);
// //     StateMatrix gradY = StateMatrix::Zero(nVar, nElem);


// //     // GRADIENT CONSTRUCTION =====================
// //     // ===============================
// //     // Interior Faces
// //     // ===============================
// //     for (int f_i = 0; f_i < I2E.rows(); ++f_i)
// //     {
// //         // Get left and right cell indices for the edge
// //         int elemL = I2E(f_i,0)-1;
// //         int elemR = I2E(f_i,2)-1;
// //         int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
// //         int faceID = mesh.elem(elemL)._faceID[lfL];  // global face index

// //         // Get edge normal vector and length
// //         Eigen::Vector2d normal = In.row(f_i);
// //         double edge_length = mesh.length(faceID);

// //         // Get cell center state values for left and right cells
// //         /** Instead of storing a copy of the state vector, we can directly reference the elements */
// //         const State& UL = U.col(elemL);
// //         const State& UR = U.col(elemR);

// //         // Compute state value at edge midpoint as average of left and right cell center values
// //         State Uhat = 0.5 * (UL + UR);

// //         // Compute gradient contribution from edge and add to left cell gradient and subtract from right cell gradient
// //         gradX.col(elemL) += Uhat * normal(0) * edge_length;
// //         gradY.col(elemL) += Uhat * normal(1) * edge_length;
// //         gradX.col(elemR) -= Uhat * normal(0) * edge_length; // opposite sign for right cell
// //         gradY.col(elemR) -= Uhat * normal(1) * edge_length; // opposite sign for right cell
// //     }

// //     // ===============================
// //     // Boundary Faces
// //     // ===============================
// //     for (int f_b = 0; f_b < B2E.rows(); ++f_b)
// //     {
// //         // Get left cell index for the edge (there is no right cell since this is a boundary edge)
// //         int elem   = B2E(f_b,0)-1;
// //         int lf     = B2E(f_b,1) - 1;
// //         int faceID = mesh.elem(elem)._faceID[lf];
// //         int bgroup = B2E(f_b,2);

// //         // Get edge normal vector and length
// //         Eigen::Vector2d normal = Bn.row(f_b);
// //         double edge_length = mesh.length(faceID);

// //         // Get cell center state values for left cell
// //         const State& UL = U.col(elem);

// //         // Use boundary state based on bgroup
// //         State Ub = State::Zero(); // placeholder, need to compute boundary state based on bgroup and interior cell center state
// //         // bunch of if statements here to compute the boundary state based on the bgroup and the interior cell center state

// //         // Compute gradient contribution from edge and add to left cell gradient
// //         gradX.col(elem) += Ub * normal(0) * edge_length;
// //         gradY.col(elem) += Ub * normal(1) * edge_length;
// //     }

// //     // ===============================
// //     // Divide gradient by Area
// //     // ===============================
// //     for (int e = 0; e < nElem; ++e)
// //     {
// //         gradX.col(e) /= Area(e);
// //         gradY.col(e) /= Area(e);
// //     }

// //     // // END OF GRADIENT CONSTRUCTION =====================


// //     // Implement limiting

// //     // RESIDUAL CONSTRUCTION =====================
// //     // ===============================
// //     // Interior Faces
// //     // ===============================
// //     for (int f_i = 0; f_i < I2E.rows(); ++f_i)
// //     {
// //         // Get left and right cell indices for the edge
// //         int elemL = I2E(f_i,0)-1;
// //         int elemR = I2E(f_i,2)-1;
// //         int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
// //         int faceID = mesh.elem(elemL)._faceID[lfL];  // global face index

// //         // Get edge normal vector, length and midpoint coordinates
// //         Eigen::Vector2d normal = In.row(f_i);
// //         double edge_length = mesh.length(faceID);
// //         const auto& face = mesh.face(faceID);
// //         Eigen::Vector2d midpoint = 0.5 * (mesh.node(face._pointID[0]) + mesh.node(face._pointID[1]));

// //         // Get state vector at the edge midpoint for left face
// //         const State& UL = U.col(elemL);
// //         Eigen::Vector2d centroidL = mesh.centroid(elemL);
// //         Eigen::Vector2d dxL = midpoint - centroidL;
// //         State UhatL = UL + gradX.col(elemL) * dxL(0) + gradY.col(elemL) * dxL(1);

// //         // Get state vector at the edge midpoint for right face
// //         const State& UR = U.col(elemR);
// //         Eigen::Vector2d centroidR = mesh.centroid(elemR);
// //         Eigen::Vector2d dxR = midpoint - centroidR;
// //         State UhatR = UR + gradX.col(elemR) * dxR(0) + gradY.col(elemR) * dxR(1);

// //         // use flux function to get the flux at the edge midpoint using uL and uR
// //         Eigen::Vector4d flux = numFlux_(UhatL, UhatR, gamma_, normal);
// //         flux *= edge_length;

// //         // add the flux contribution to the residual of the left cell and subtract it from the residual of the right cell. 
// //             // When this is done over all the edges this will construct the full residual over all the cells in the domain
// //         R.col(elemL) += flux;
// //         R.col(elemR) -= flux;   // opposite sign

// //         // compute the wave speed at the edge and update the maximum wave speed for the left and right cells if necessary
// //     }

// //     // ===============================
// //     // Boundary Faces
// //     // ===============================
// //     for (int f_b = 0; f_b < B2E.rows(); ++f_b)
// //     {
// //         // Get left cell index for the edge (there is no right cell since this is a boundary edge)
// //         int elem   = B2E(f_b,0)-1;
// //         int lf     = B2E(f_b,1) - 1;
// //         int faceID = mesh.elem(elem)._faceID[lf];
// //         int bgroup = B2E(f_b,2);

// //         // Get edge normal vector and length
// //         Eigen::Vector2d normal = Bn.row(f_b);
// //         double edge_length = mesh.length(faceID);

// //         // REST NEEDS TO BE DONE
// //     }

// //     return R;
// // }

// // } // namespace solver