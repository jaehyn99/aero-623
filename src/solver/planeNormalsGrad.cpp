// Implementing only the gradient function plane normals reconstruction

# include "solver/planeNormalsGrad.h"
# include "mesh/TriangularMesh.h"

# include <Eigen/Dense>

using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using State = Eigen::Vector4d;

namespace solver {
// ======================================================
// GRADIENT CONSTRUCTION FUNCTION
// ======================================================

void SecondOrderEuler::computeGrad_pn(
    TriangularMesh& mesh,
    const Eigen::MatrixXi& I2E,
    const Eigen::MatrixXi& B2E,
    const Eigen::MatrixXd& In,
    const Eigen::MatrixXd& Bn,
    const Eigen::VectorXd& Area,
    const StateMatrix& U,
    StateMatrix& gradX,
    StateMatrix& gradY) const
{
    // ===============================
    // Interior Faces
    // ===============================
    for (int f_i = 0; f_i<I2E.rows(); ++f_i)
    {
        // Get left and right cell indices for the edge
        int elemL = I2E(f_i,0)-1;
        int elemR = I2E(f_i,2)-1;
        int lfL    = I2E(f_i,1) - 1;                 // local face on elemL
        int faceID = mesh.elem(elemL)._faceID[lfL];  // global face index

        // Get edge normal vector and length
        Eigen::Vector2d normal = In.row(f_i);
        double edge_length = mesh.length(faceID);

        // Get cell center state values for left and right cells
        /** Instead of storing a copy of the state vector, we can directly reference the elements */
        const State& UL = U.col(elemL);
        const State& UR = U.col(elemR);

        // Compute state value at edge midpoint as average of left and right cell center values
        State Uhat = 0.5 * (UL + UR);

        // Compute gradient contribution from edge and add to left cell gradient and subtract from right cell gradient
        gradX.col(elemL) += Uhat * normal(0) * edge_length;
        gradY.col(elemL) += Uhat * normal(1) * edge_length;
        gradX.col(elemR) -= Uhat * normal(0) * edge_length; // opposite sign for right cell
        gradY.col(elemR) -= Uhat * normal(1) * edge_length; // opposite sign for right cell
    }

    // ===============================
    // Boundary Faces
    // ===============================
    for (int f_b = 0; f_b<B2E.rows(); ++f_b)
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

        // Construct p vectors for plane normal gradient construction
        // What are the neighnoring cells for this boundary cell?
        for (int k = 0; k<4; ++k)
        {
            // Eigen::Vector2d xe = centroid[e];
            // Eigen::Vector2d x1 = centroid[n1];
            // Eigen::Vector2d x2 = centroid[n2];

            // double u0 = U[e](k);
            // double u1 = U[n1](k);
            // double u2 = U[n2](k);

            // Eigen::Vector3d p0(xe(0), xe(1), u0);
            // Eigen::Vector3d p1(x1(0), x1(1), u1);
            // Eigen::Vector3d p2(x2(0), x2(1), u2);
            
            // // Compute plane normal using cross product of (p1-p0) and (p2-p0)
            // Eigen::Vector3d plane_normal = (p1 - p0).cross(p2 - p0);

            // NOT SURE HOW TO COMPUTE THE GRADIENT NOW
        }


    }
}

} //namespace solver