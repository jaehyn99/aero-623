// Implementation of the gradient function using walk over 3 edges for interior cells and plane normals reconstruction for boundaries

# include "hybridWalkPNGrad.h"
# include "TriangularMesh.h"
# include "StateMesh.h"

using State = Eigen::Vector4d;

// ======================================================
// GRADIENT CONSTRUCTION FUNCTION
// ======================================================

std::vector<Eigen::Matrix<double,4,2>> HybridWalkPNGrad::computeGradient(const StateMesh& U) const
{
    Eigen::MatrixXd gradX = Eigen::MatrixXd::Zero(U.stateCount(), U.cellCount());
    Eigen::MatrixXd gradY = Eigen::MatrixXd::Zero(U.stateCount(), U.cellCount());
    auto mesh = U.mesh();
    std::vector<bool> visited(mesh->numFaces(), false);

    // Interior Faces
    for (std::size_t i = 0; i < mesh->numFaces(); i++){
        if (visited[i]) continue; // Skip if this face has already been visited (handles periodic pairs)
        const auto& face = mesh->face(i);
        visited[i] = true;
        if (face.isPeriodicFace())
            visited[face._periodicFaceID] = true; // Mark the periodic counterpart as visited to avoid double counting

        if (!face.isBoundaryFace()){
            // Get left and right cell indices for the edge
            int elemL  = face._elemID[0];
            int elemR  = face.isPeriodicFace() ? face._periodicElemID : face._elemID[1];

            // Get edge normal vector and length
            Eigen::Vector2d normal = face._normal;
            double edge_length = mesh->length(i);
            // if (face.isPeriodicFace()){
            //     if (face._elemID[0] > face._periodicElemID) normal *= -1;
            // }

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
        }
    }

    // ===============================
    // Divide gradient by Area
    // ===============================
    for (int e = 0; e < U.cellCount(); ++e)
    {
        gradX.col(e) /= mesh->area(e);
        gradY.col(e) /= mesh->area(e);
    }

    // ===============================
    // Boundary Faces
    // ===============================
    for (std::size_t b = 0; b < mesh->numFaces(); b++){
        const auto& face = mesh->face(b);
        if (!face.isBoundaryFace()) continue;

        // Get left cell index for the edge (there is no right cell since this is a boundary edge)
        int elem   = face._elemID[0];

        // // Get edge normal vector and length
        // Eigen::Vector2d normal = face._normal;
        // double edge_length = mesh->length(b);

        // // Get cell center state values for left cell
        // const State& UL = U.cell(elem);

        std::array<int,2> neighbors;
        int count = 0;
        const auto& elemObj = mesh->elem(elem);
        for (std::size_t j = 0; j < 3; ++j){
            int neighbor_fID = elemObj._faceID[j];

            // Skip the boundary face we are currently on
            if (neighbor_fID == b) continue;
            if (neighbor_fID == face._periodicFaceID) continue;

            const auto& f = mesh->face(neighbor_fID);
            int n0 = f._elemID[0];
            int n1 = f._elemID[1];
            int neigh = -1;

            if (n0 == elem) neigh = n1;
            else if (n1 == elem) neigh = n0;

            // Only store valid interior neighbor
            if (neigh != -1)
                neighbors[count++] = neigh;
        }

        // Construct p vectors for plane normal gradient construction
        // What are the neighnoring cells for this boundary cell?
        Eigen::Vector2d xe = mesh->centroid(elem);
        Eigen::Vector2d x1 = mesh->centroid(neighbors[0]);
        Eigen::Vector2d x2 = mesh->centroid(neighbors[1]);

        for (int k = 0; k<4; ++k)
        {
            double u0 = U(k, elem);
            double u1 = U(k, neighbors[0]);
            double u2 = U(k, neighbors[1]);

            Eigen::Vector3d p0(xe(0), xe(1), u0);
            Eigen::Vector3d p1(x1(0), x1(1), u1);
            Eigen::Vector3d p2(x2(0), x2(1), u2);

            // Compute plane normal using cross product of (p1-p0) and (p2-p0)
            Eigen::Vector3d plane_normal = (p1 - p0).cross(p2 - p0);

            // Check if plane normal has 0 in 3rd component to avoid division by zero
            if (std::abs(plane_normal(2)) < 1e-15) plane_normal(2) = 1e-15;
            gradX(k, elem) = -plane_normal(0) / plane_normal(2);
            gradY(k, elem) = -plane_normal(1) / plane_normal(2);
        }
    }

    std::vector<Eigen::Matrix<double,4,2>> Lgrad(U.cellCount());
    for (int e = 0; e < U.cellCount(); ++e)
    {
        Lgrad[e].col(0) = gradX.col(e);   // dx
        Lgrad[e].col(1) = gradY.col(e);   // dy
    }

    return Lgrad;
}

// // Implementation of the gradient function using walk over 3 edges for interior cells and plane normals reconstruction for boundaries

// # include "HybridWalkPNGrad.h"
// # include "TriangularMesh.h"
// # include "StateMesh.h"

// // not sure if  i should keep this ======================
// #include <Eigen/StdVector>
// // ======================================================

// # include <Eigen/Dense>

// // ======================================================
// // Dubug Includes
// #include <fstream>
// #include <iomanip>
// #include <cstdlib>
// // ======================================================

// using State = Eigen::Vector4d;

// // ======================================================
// // GRADIENT CONSTRUCTION FUNCTION
// // ======================================================

// std::vector<Eigen::Matrix<double,4,2>> HybridWalkPNGrad::computeGradient(const Eigen::MatrixXi& I2E,
//                                                                         const Eigen::MatrixXi& B2E,
//                                                                         const Eigen::MatrixXd& In,
//                                                                         const Eigen::MatrixXd& Bn,
//                                                                         const Eigen::VectorXd& Area,
//                                                                         const StateMesh& U) const
// {

//     // open once per call (overwrite each run)
//     std::ofstream dbg("periodic_edges_debug.txt", std::ios::trunc);
//     dbg << std::setprecision(17);
//     dbg << "Periodic edge debug log\n";

//     Eigen::MatrixXd gradX = Eigen::MatrixXd::Zero(U.stateCount(), U.cellCount());
//     Eigen::MatrixXd gradY = Eigen::MatrixXd::Zero(U.stateCount(), U.cellCount());

//     auto mesh = U.mesh();

//     std::vector<bool> visited(mesh->numFaces(), false);


//     // ===============================
//     // Interior Faces
//     // ===============================
//     for (int f_i = 0; f_i<I2E.rows(); ++f_i)
//     {
//         // Get left and right cell indices for the edge
//         int elemL = I2E(f_i,0)-1;
//         int elemR = I2E(f_i,2)-1;
//         int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
//         int faceID = mesh->elem(elemL)._faceID[lfL];  // global face index
        
//         const auto& face = mesh->face(faceID);
//         if (visited[faceID]) {
//             continue; // Skip if this face has already been visited (handles periodic pairs)
//         }
//         visited[faceID] = true;
//         if (face.isPeriodicFace())
//             visited[face._periodicFaceID] = true; // Mark the periodic counterpart as visited to avoid double counting

//         // Get edge normal vector and length
//         Eigen::Vector2d normal = In.row(f_i);
//         double edge_length = mesh->length(faceID);

//         if (face.isPeriodicFace()){
//             if (face._elemID[0] > face._periodicElemID){
//                 normal *= -1;
//             }
//         }

//         // Get cell center state values for left and right cells
//         /** Instead of storing a copy of the state vector, we can directly reference the elements */
//         State UL = U.cell(elemL);
//         State UR = U.cell(elemR);

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
//             const Eigen::Vector2d cL = mesh->centroid(elemL);
//             const Eigen::Vector2d cR = mesh->centroid(elemR);

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
//     // Divide gradient by Area
//     // ===============================
//     for (int e = 0; e < U.cellCount(); ++e)
//     {
//         gradX.col(e) /= Area(e);
//         gradY.col(e) /= Area(e);
//     }

//     // ===============================
//     // Boundary Faces
//     // ===============================
//     for (int f_b = 0; f_b<B2E.rows(); ++f_b)
//     {
//         // Get left cell index for the edge (there is no right cell since this is a boundary edge)
//         int elem   = B2E(f_b,0)-1;
//         int lf     = B2E(f_b,1) - 1;
//         int faceID = mesh->elem(elem)._faceID[lf];
//         int bgroup = B2E(f_b,2);

//         // Get edge normal vector and length
//         Eigen::Vector2d normal = Bn.row(f_b);
//         double edge_length = mesh->length(faceID);

//         // Get cell center state values for left cell
//         const State& UL = U.cell(elem);

//         std::array<int,2> neighbors;
//         int count = 0;

//         const auto& elemObj = mesh->elem(elem);

//         for (int j = 0; j < 3; ++j)
//         {
//             int neighbor_fID = elemObj._faceID[j];

//             // Skip the boundary face we are currently on
//             if (neighbor_fID == faceID)
//                 continue;

//             const auto& f = mesh->face(neighbor_fID);

//             int n0 = f._elemID[0];
//             int n1 = f._elemID[1];

//             int neigh = -1;

//             if (n0 == elem) neigh = n1;
//             else if (n1 == elem) neigh = n0;

//             // Only store valid interior neighbor
//             if (neigh != -1)
//                 neighbors[count++] = neigh;
//         }

//         // Construct p vectors for plane normal gradient construction
//         // What are the neighnoring cells for this boundary cell?
//         Eigen::Vector2d xe = mesh->centroid(elem);
//         Eigen::Vector2d x1 = mesh->centroid(neighbors[0]);
//         Eigen::Vector2d x2 = mesh->centroid(neighbors[1]);

//         for (int k = 0; k<4; ++k)
//         {
//             double u0 = U(k, elem);
//             double u1 = U(k, neighbors[0]);
//             double u2 = U(k, neighbors[1]);

//             Eigen::Vector3d p0(xe(0), xe(1), u0);
//             Eigen::Vector3d p1(x1(0), x1(1), u1);
//             Eigen::Vector3d p2(x2(0), x2(1), u2);

//             // Compute plane normal using cross product of (p1-p0) and (p2-p0)
//             Eigen::Vector3d plane_normal = (p1 - p0).cross(p2 - p0);

//             // Check if plane normal has 0 in 3rd component to avoid division by zero
//             if (std::abs(plane_normal(2)) < 1e-15)
//             {
//                 plane_normal(2) = 1e-15;
//             }
//             gradX(k , elem) = - plane_normal(0) / plane_normal(2);
//             gradY(k , elem) = - plane_normal(1) / plane_normal(2);
//         }
//     }

//     std::vector<Eigen::Matrix<double,4,2>> Lgrad;
//     Lgrad.resize(U.cellCount());

//     for (int e = 0; e < U.cellCount(); ++e)
//     {
//         Lgrad[e].col(0) = gradX.col(e);   // dx
//         Lgrad[e].col(1) = gradY.col(e);   // dy
//     }

//     return Lgrad;
// }