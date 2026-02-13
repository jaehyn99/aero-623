// WHAT I HAVE IMPLEMENETED HERE IS USING PLANE NORMALS FOR GRADIENT CONSTRUCTION AT BOUNDARY CELLS

# include "solver/SecondorderEuler.h"
# include "mesh/TriangularMesh.h"

# include <Eigen/Dense>

using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using State = Eigen::Vector4d;
using Grad  = Eigen::Matrix<double,4,2>;


StateMatrix
computeResidual(const Eigen::MatrixXi& I2E,
                const Eigen::MatrixXi& B2E,
                const Eigen::MatrixXd& In,
                const Eigen::MatrixXd& Bn,
                const Eigen::VectorXd& Area,
                const StateMatrix& U)
{
    // Get number of elements and number of variables from the input state vector
    const int nElem = U.cols();
    const int nVar  = U.rows();

    // initialize the residual on all cells to zero
    StateMatrix R(nVar, nElem);

    // initialize the gradient on all cells to zero
    StateMatrix gradX = StateMatrix::Zero(nVar, nElem);
    StateMatrix gradY = StateMatrix::Zero(nVar, nElem);


    // GRADIENT CONSTRUCTION =====================
    // ===============================
    // Interior Faces
    // ===============================
    for (int f_i = 0; f_i < I2E.rows(); ++f_i)
    {
        // Get left and right cell indices for the edge
        int elemL = I2E(f_i,0)-1;
        int elemR = I2E(f_i,2)-1;

        // Get edge normal vector and length
        Eigen::Vector2d normal = In.row(f_i);
        int edge_length = 1; // placeholder, need to compute edge length from node coordinates

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
    for (int f_b = 0; f_b < B2E.rows(); ++f_b)
    {
        // Get left cell index for the edge (there is no right cell since this is a boundary edge)
        int elem   = B2E(f_b,0)-1;
        int bgroup = B2E(f_b,2);

        // Get edge normal vector and length
        Eigen::Vector2d normal = Bn.row(f_b);
        int edge_length = 1; // placeholder, need to compute edge length from node coordinates

        // Get cell center state values for left cell
        const State& UL = U.col(elem);

        // Construct p vectors for plane normal gradient construction
        // needs an E2E (element to element) connectivity to get the neighboring cell for the boundary cells
        // construct p vector as [x, y, u].T where x, y are coords of cell center and u is each element of the state vector at cell center.
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

    // ===============================
    // Divide gradient by Area
    // ===============================
    for (int e = 0; e < nElem; ++e)
    {
        gradX.col(e) /= Area(e);
        gradY.col(e) /= Area(e);
    }

    // END OF GRADIENT CONSTRUCTION =====================


    // RESIDUAL CONSTRUCTION =====================
    // ===============================
    // Interior Faces
    // ===============================
    for (int f_i = 0; f_i < I2E.rows(); ++f_i)
    {
        // Get left and right cell indices for the edge
        int elemL = I2E(f_i,0);
        int elemR = I2E(f_i,2);

        // Get edge normal vector, length and midpoint coordinates
        Eigen::Vector2d normal = In.row(f_i);
        int edge_length = 1; // placeholder, need to compute edge length from node coordinates
        Eigen::Vector2d midpoint; // placeholder, need to compute midpoint from node coordinates

        // Get state vector at the edge midpoint for left face
        const State& UL = U[elemL];
        const Grad& gradL = grad[elemL];
        State uL = UL + gradL * (midpoint - Eigen::Vector2d::Zero()); // placeholder for cell center coordinates

        // Get state vector at the edge midpoint for right face
        const State& UR = U[elemR];
        const Grad& gradR = grad[elemR];
        State uR = UR + gradR * (midpoint - Eigen::Vector2d::Zero()); // placeholder for cell center coordinates

        // use flux function to get the flux at the edge midpoint using uL and uR
        State flux = computeNumericalFlux(uL, uR, normal);

        // add the flux contribution to the residual of the left cell and subtract it from the residual of the right cell. 
            // When this is done over all the edges this will construct the full residual over all the cells in the domain
        R[elemL] += flux;
        R[elemR] -= flux;   // opposite sign

        // compute the wave speed at the edge and update the maximum wave speed for the left and right cells if necessary
    }

    // ===============================
    // Boundary Faces
    // ===============================
    for (int f_b = 0; f_b < B2E.rows(); ++f_b)
    {
        // Get left cell index for the edge (there is no right cell since this is a boundary edge)
        int elem   = B2E(f_b,0);
        int bgroup = B2E(f_b,2);

        // Get edge normal vector and length
        Eigen::Vector2d normal = Bn.row(f_b);
        int edge_length = 1; // placeholder, need to compute edge length from node coordinates

        // REST NEEDS TO BE DONE
    }

    return R;
}