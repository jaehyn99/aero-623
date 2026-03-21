// # include "solver/FirstOrderEuler_alt.h"
// # include "mesh/TriangularMesh.h"
// #include <iostream>
// #include <cmath>


// # include <Eigen/Dense>

// using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
// using State = Eigen::Vector4d;

// namespace solver {
// // ======================================================
// // RESIDUAL CONSTRUCTION FUNCTION
// // ======================================================
// ResidualResult FirstOrderEuler_alt::computeResidual(
//     TriangularMesh& mesh,
//     const Eigen::MatrixXi& I2E,
//     const Eigen::MatrixXi& B2E,
//     const Eigen::MatrixXd& In,
//     const Eigen::MatrixXd& Bn,
//     const Eigen::VectorXd& Area,
//     const StateMatrix& U) const
// {
//     std::vector<bool> visited(mesh.numFaces(), false);

//     StateMatrix R = StateMatrix::Zero(U.rows(), U.cols());
//     Eigen::VectorXd waveSum   = Eigen::VectorXd::Zero(U.cols());
//     Eigen::VectorXd perimeter   = Eigen::VectorXd::Zero(U.cols());

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
//         if (face._periodicFaceID != -1)
//         {
//             visited[face._periodicFaceID] = true; // Mark the periodic counterpart as visited to avoid double counting
//         }

//         // Get edge normal vector, length and midpoint coordinates
//         Eigen::Vector2d normal = In.row(f_i);
//         double edge_length = mesh.length(faceID);
//         if (face.isPeriodicFace()){
//             if (face._elemID[0] > face._periodicElemID){
//                 normal *= -1;
//             }
//         }

//         // Get state vector at the edge midpoint for left face
//         int faceIDL = mesh.elem(elemL)._faceID[lfL];
//         const auto& faceL = mesh.face(faceIDL);
//         const State& UL = U.col(elemL);


//         // Get state vector at the edge midpoint for right face
//         int lfR = I2E(f_i,3) - 1; // local face on elemR
//         int faceIDR = mesh.elem(elemR)._faceID[lfR];
//         const auto& faceR = mesh.face(faceIDR);
//         const State& UR = U.col(elemR);

//         // use flux function to get the flux at the edge midpoint using uL and uR
//         auto result = numFlux_(UL, UR, gamma_, normal);
//         Eigen::Vector4d flux = result.flux;
//         double lambda_max    = result.maxLambda;
//         flux *= edge_length;

//         if (!std::isfinite(lambda_max))
//         {
//             std::cout << "\nNaN lambda at face " << faceID << "\n";
//             std::cout << "elemL: " << elemL << " elemR: " << elemR << "\n";
//             std::cout << "UL: " << UL.transpose() << "\n";
//             std::cout << "UR: " << UR.transpose() << "\n";
//             std::cout << "normal: " << normal.transpose() << "\n";
//             std::cout << "edge_length: " << edge_length << "\n";

//             throw std::runtime_error("NaN lambda from flux");
//         }

//         // add the flux contribution to the residual of the left cell and subtract it from the residual of the right cell. 
//             // When this is done over all the edges this will construct the full residual over all the cells in the domain
//         R.col(elemL) += flux;
//         R.col(elemR) -= flux;   // opposite sign

//         // compute the wave speed at the edge and update the maximum wave speed for the left and right cells if necessary
//         double contrib = lambda_max * edge_length;
//         waveSum(elemL) += contrib;
//         waveSum(elemR) += contrib;

//         if (!std::isfinite(waveSum(elemL)) || !std::isfinite(waveSum(elemR)))
//         {
//             std::cout << "\nNaN lambda at face " << faceID << "\n";
//             std::cout << "elemL: " << elemL << " elemR: " << elemR << "\n";
//             std::cout << "UL: " << UL.transpose() << "\n";
//             std::cout << "UR: " << UR.transpose() << "\n";
//             std::cout << "normal: " << normal.transpose() << "\n";
//             std::cout << "edge_length: " << edge_length << "\n";

//             throw std::runtime_error("NaN lambda from flux");
//         }

//         perimeter(elemL) += edge_length;
//         perimeter(elemR) += edge_length;
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
//             auto result = inletFlux_(UL, gamma_, normal);
//             Eigen::Vector4d flux = result.flux;
//             double lambda_max = result.maxLambda;
//             flux *= edge_length;
//             R.col(elem) += flux; // add contribution to residual of left cell
//             waveSum(elem) += lambda_max * edge_length; // add contribution to wave speed for left cell
//             perimeter(elem) += edge_length;

//             // DEBUG BLOCK
//             if (!result.flux.allFinite() || !std::isfinite(result.maxLambda))
//             {
//                 std::cout << "\n--- Boundary Debug ---\n";
//                 std::cout << "Face: " << faceID
//                         << "  Cell: " << elem
//                         << "  bgroup: " << bgroup << "\n";

//                 std::cout << "Flux components:\n";

//                 for (int k = 0; k < 4; ++k)
//                 {
//                     if (std::isnan(flux(k)))
//                         std::cout << "  F(" << k << ") = NaN\n";
//                     else if (!std::isfinite(flux(k)))
//                         std::cout << "  F(" << k << ") = Inf\n";
//                     else
//                         std::cout << "  F(" << k << ") = " << flux(k) << "\n";
//                 }

//                 if (std::isnan(result.maxLambda))
//                     std::cout << "Lambda = NaN\n";
//                 else if (!std::isfinite(result.maxLambda))
//                     std::cout << "Lambda = Inf\n";
//                 else
//                     std::cout << "Lambda = " << result.maxLambda << "\n";
//             }
//         }
//         else if (bgroup == 7) { // outlet
//             auto result = outletFlux_(UL, gamma_, normal);
//             Eigen::Vector4d flux = result.flux;
//             double lambda_max = result.maxLambda;
//             flux *= edge_length;
//             R.col(elem) += flux; // add contribution to residual of left cell
//             waveSum(elem) += lambda_max * edge_length; // add contribution to wave speed for left cell
//             perimeter(elem) += edge_length;

//             // DEGUG BLOCK
//             if (!result.flux.allFinite() || !std::isfinite(result.maxLambda))
//             {
//                 std::cout << "NaN from boundary face " << faceID
//                         << " in cell " << elem
//                         << " bgroup: " << bgroup << "\n";
//             }
//         }
//         else if (bgroup == 1 || bgroup == 5) { // wall
//             auto result = wallFlux_(UL, gamma_, normal);
//             Eigen::Vector4d flux = result.flux;
//             double lambda_max = result.maxLambda;
//             flux *= edge_length;
//             R.col(elem) += flux; // add contribution to residual of left cell
//             waveSum(elem) += lambda_max * edge_length; // add contribution to wave speed for left cell
//             perimeter(elem) += edge_length;

//             // DEGUG BLOCK
//             if (!result.flux.allFinite() || !std::isfinite(result.maxLambda))
//             {
//                 std::cout << "NaN from boundary face " << faceID
//                         << " in cell " << elem
//                         << " bgroup: " << bgroup << "\n";
//             }
//         }
//         else
//             throw std::runtime_error("Unknown boundary group in computeResidual");

//     }

//     Eigen::VectorXd waveSpeed(U.cols());

//     for (int i = 0; i < U.cols(); ++i)
//     {
//         waveSpeed(i) = waveSum(i) / perimeter(i);
//     }

//     return {R, waveSpeed, perimeter};
// }

// // ======================================================
// // Time Step Function
// // ======================================================
// void FirstOrderEuler_alt::marchToSteadyState(
//     TriangularMesh& mesh,
//     const Eigen::MatrixXi& I2E,
//     const Eigen::MatrixXi& B2E,
//     const Eigen::MatrixXd& In,
//     const Eigen::MatrixXd& Bn,
//     const Eigen::VectorXd& Area,
//     StateMatrix& U,
//     double CFL,
//     int maxIter,
//     double tol) const
// {
//     for (int iter = 0; iter < maxIter; ++iter)
//     {
//         // --------------------------------------
//         // Compute residual and wave speeds
//         // --------------------------------------
//         auto result = computeResidual(mesh, I2E, B2E, In, Bn, Area, U);

//         StateMatrix R = result.R;
//         Eigen::VectorXd lambda = result.waveSpeed;
//         Eigen::VectorXd perimeter = result.perimeter;

//         // --------------------------------------
//         // Compute global time step
//         // --------------------------------------
//         double dt = 1e2;
//         Eigen::VectorXd dt_vec = Eigen::VectorXd::Zero(U.cols());

//         for (int i = 0; i < U.cols(); ++i)
//         {

//             // DEBUG BLOCK
//             if (!std::isfinite(lambda(i)))
//             {
//                 std::cout << "\n==== DEBUG CELL " << i << " ====\n";

//                 std::cout << "Area: " << Area(i) << "\n";
//                 std::cout << "Perimeter: " << perimeter(i) << "\n";
//                 std::cout << "waveSum: " << result.waveSpeed(i) * perimeter(i) << "\n";

//                 Eigen::Vector2d c = mesh.centroid(i);
//                 std::cout << "Centroid: " << c.transpose() << "\n";

//                 // Print faces of this element
//                 std::cout << "Faces:\n";

//                 for (int lf = 0; lf < 3; ++lf)
//                 {
//                     int faceID = mesh.elem(i)._faceID[lf];
//                     const auto& face = mesh.face(faceID);

//                     std::cout << "  Local face " << lf
//                             << "  Global faceID " << faceID
//                             << "  length: " << mesh.length(faceID)
//                             << "  periodic: " << face._periodicFaceID
//                             << "\n";

//                     // -------------------------------------------------
//                     // Determine if this face is interior or boundary
//                     // -------------------------------------------------

//                     bool foundInterior = false;
//                     int neighbor = -1;

//                     for (int f_i = 0; f_i < I2E.rows(); ++f_i)
//                     {
//                         int elemL = I2E(f_i,0) - 1;
//                         int elemR = I2E(f_i,2) - 1;

//                         if (elemL == i && mesh.elem(elemL)._faceID[I2E(f_i,1)-1] == faceID)
//                         {
//                             neighbor = elemR;
//                             foundInterior = true;
//                             break;
//                         }
//                         if (elemR == i && mesh.elem(elemR)._faceID[I2E(f_i,3)-1] == faceID)
//                         {
//                             neighbor = elemL;
//                             foundInterior = true;
//                             break;
//                         }
//                     }

//                     if (foundInterior)
//                     {
//                         std::cout << "     Interior face\n";
//                         std::cout << "     Neighbor cell: " << neighbor << "\n";

//                         Eigen::Vector2d nc = mesh.centroid(neighbor);
//                         std::cout << "     Neighbor centroid: " << nc.transpose() << "\n";

//                         std::cout << "     U(i): " << U.col(i).transpose() << "\n";
//                         std::cout << "     U(neighbor): " << U.col(neighbor).transpose() << "\n";
//                     }
//                     else
//                     {
//                         std::cout << "     Boundary face\n";
//                     }
//                 }

//             }
//             // --------------------------------------

            
//             double d_i = 2 * Area(i) / perimeter(i);
//             dt_vec(i) = CFL * d_i / lambda(i);
//         }
//         double dt_glob = dt_vec.minCoeff();
//         dt = std::min(dt, dt_glob);
//         // --------------------------------------
//         // Compute normalized L2 residual
//         // --------------------------------------
//         double resNorm = 0.0;

//         for (int i = 0; i < U.cols(); ++i)
//         {
//             for (int v = 0; v < 4; ++v)
//             {
//                 double r = R(v,i);
//                 resNorm += r*r;
//             }
//         }

//         resNorm = std::sqrt(resNorm);

//         // Normalize by number of cells
//         resNorm /= U.cols();

//         std::cout << "Iter: " << iter
//                   << "  Residual L2: " << resNorm
//                   << "  log10: " << std::log10(resNorm)
//                   << "  dt: " << dt
//                   << std::endl;

//         // --------------------------------------
//         // Convergence check
//         // --------------------------------------
//         if (resNorm < tol)
//         {
//             std::cout << "Converged." << std::endl;
//             break;
//         }

//         // --------------------------------------
//         // Forward Euler update
//         // --------------------------------------
//         for (int i = 0; i < U.cols(); ++i)
//         {
//             U.col(i) -= dt * R.col(i) / Area(i);
//         }
//     }
// }


// } // namespace solver