
// Get edge normals and cell areas for second order euler solver
  
// Get cell center values from first order solver

// start do while loop here to reduce residual below a threshold - this is steady state version

    // initialize gradient to zero on all cells

    // loop over all edges in the entire flow domain and compute gradient contribution from each edge
        // THIS IS FOR INTERNAL EDGES
        // get left and right cell indices for the edge
        // get edge normal vector and length
        // get cell center values for left and right cells
        // Compute the u_hat value at the edge as the average of the left and right cell center values
        // Compute gradient contribution as u_hat * normal vector * edge length and add this to the gradient of the left cell
        // Subtract the same contribution from the gradient of the right cell

        // THIS IS FOR BOUNDARY EDGES
    // end for loop over edges

    // Divide each cell's gradient by its area to get the final gradient value for each cell

    // Set (initialize) the residual to zero and set the wave speed on each cell to zero

    // loop over all edges in the entire flow domain and compute the flux contribution from each edge
        // THIS IS FOR INTERNAL EDGES
        // get left and right cell indices for the edge
        // get edge normal vector, length and midpoint coordinates
        // get state vector at the edge midpoint for left face using uL(x_mid) = u(left cell center) + grad(u(left cell center)) dot (x_mid - x(left cell center))
        // uR(x_mid) = u(right cell center) + grad(u(right cell center)) dot (x_mid - x(right cell center))
        // use flux function to get the flux at the edge midpoint using uL and uR
        // add the flux contribution to the residual of the left cell and subtract it from the residual of the right cell
        // compute the wave speed at the edge and update the maximum wave speed for the left and right cells if necessary

        // THIS IS FOR BOUNDARY EDGES

    // end for loop over edges

    // Compute time step using CFL condition based on the maximum wave speed in each cell and the cell area

    // compute next state using time integration








// This is actual second order flux residual function

// initialize gradient to zero on all cells
/** Gradient is defined over all cells */

    // loop over all edges in the entire flow domain and compute gradient contribution from each edge to the left and right cells
        // THIS IS FOR INTERNAL EDGES loop over I2E for internal edges
        // get left and right cell indices for the edge -

        // get edge normal vector and length -////
        /** The edge normals and lengths are defined over all edges */

        // get cell center values for left and right cells -
        /** The cell center values are defined over all cells */

        // Compute the u_hat value at the edge as the average of the left and right cell center values -
        /** The u_hat values are defined over all edges */

        // Compute gradient contribution as u_hat * normal vector * edge length and add this to the gradient of the left cell -

        // Subtract the same contribution from the gradient of the right cell -

        // THIS IS FOR BOUNDARY EDGES loop over B2E for boundary edges
        // get left cell index for the edge (there is no right cell since this is a boundary edge)
        // get edge normal vector and length
        // get cell center values for left cell
        // use boundary state
    // end for loop over edges

    // Divide each cell's gradient by its area to get the final gradient value for each cell
    /** The cell areas are defined over all cells */

    // Set (initialize) the residual to zero and set the wave speed on each cell to zero

    // loop over  all edges in the entire flow domain and compute the flux contribution from each edge
        // THIS IS FOR INTERNAL EDGES
        // get left and right cell indices for the edge

        // get edge normal vector, length and midpoint coordinates
        /** The edge normals, lengths and midpoints are defined over all edges */

        // get state vector at the edge midpoint for left face using uL(x_mid) = u(left cell center) + grad(u(left cell center)) dot (x_mid - x(left cell center))
        
        // uR(x_mid) = u(right cell center) + grad(u(right cell center)) dot (x_mid - x(right cell center))
        /** uL and uR are defined at each edge */

        // use flux function to get the flux at the edge midpoint using uL and uR
        

        // THIS IS FOR BOUNDARY EDGES (normal vector always points outwards from the domain)
        // if inflow boundary
            // Compute exterior static temp, static P, static density, speed of sound, exterior velocity and exterior total energy using BCs
            // BC state serves as uR and the interior cell center value serves as uL
            // use flux function to get the flux at the edge midpoint using uL and uR
            // (not sure if this is correct or if F_b is computed using the interior state via Reimann Invariant and the Tt, Pt and alpha from the BCs)
        // if outflow boundary
            // use 3.3.7 to compute exterior density, 3.3.8 for b normal velocity, 3.3.9 for b tangential velocity and b total energy using b Pressure
            // BC state serves as uR and the interior cell center value serves as uL
            // use flux function to get the flux at the edge midpoint using uL and uR
            // (not sure if this is correct or if F_b is computed using the interior state via Reimann Invariant and the specified b Pressure)
        // if wall boundary
            // b Pressure is computed using 3.3.4, where tangential b velocity is computed using 3.3.5
            // Boundary flux is computed using this b Pressure as F_b = [0, P_b * n_x, P_b * n_y, 0] where n_x and n_y are components of normal vector

        // RESIDUAL CONSTRUCTION
        // add the flux contribution to the residual of the left cell and subtract it from the residual of the right cell. 
            // When this is done over all the edges this will construct the full residual over all the cells in the domain
        /** The residuals are defined over all cells */
        // compute the wave speed at the edge and update the maximum wave speed for the left and right cells if necessary
    // end for loop over edges


// Have to loop over I2E for internal edges and B2E for boundary edges to compute the flux contributions and construct the residuals

// KEEP TRACK OF THE FACT THAT INDEXINF STARTS FROM 1 IN MESH CONNECTIVITY FILES -
// if so, subtract 1 from element indices obtained from I2E and B2E when referencing the state vector defined over all cells

// WHAT I HAVE IMPLEMENTED HERE IS USING BOUNDARY STATE FOR GRADIENT CONSTRUCTION AT BOUNDARY CELLS

# include "solver/SecondorderEuler.h"
# include "mesh/TriangularMesh.h"

# include <Eigen/Dense>

using StateMatrix = Eigen::Matrix<double,4,Eigen::Dynamic>;
using State = Eigen::Vector4d;
using Grad  = Eigen::Matrix<double,4,2>;



namespace solver {

StateMatrix
SecondOrderEuler::computeResidual(const Eigen::MatrixXi& I2E,
                                  const Eigen::MatrixXi& B2E,
                                  const Eigen::MatrixXd& In,
                                  const Eigen::MatrixXd& Bn,
                                  const Eigen::VectorXd& Area,
                                  const StateMatrix& U) const
{

    // Can be replaced by other meshes
    ::TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");

    // Get number of elements and number of variables from the input state vector
    const int nElem = U.cols();
    const int nVar  = U.rows();

    // initialize the residual on all cells to zero
    StateMatrix R = StateMatrix::Zero(nVar, nElem);

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

        // Use boundary state based on bgroup
        State Ub = State::Zero(); // placeholder, need to compute boundary state based on bgroup and interior cell center state
        // bunch of if statements here to compute the boundary state based on the bgroup and the interior cell center state

        // Compute gradient contribution from edge and add to left cell gradient
        gradX.col(elem) += Ub * normal(0) * edge_length;
        gradY.col(elem) += Ub * normal(1) * edge_length;
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


    // Implement limiting

    // RESIDUAL CONSTRUCTION =====================
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

        // Get edge normal vector, length and midpoint coordinates
        Eigen::Vector2d normal = In.row(f_i);
        double edge_length = mesh.length(faceID);
        Eigen::Vector2d midpoint; // placeholder, need to compute midpoint from node coordinates

        // Get state vector at the edge midpoint for left face
        const State& UL = U.col(elemL);
        Eigen::Vector2d centroidL = mesh.centroid(elemL);
        Eigen::Vector2d dxL = midpoint - centroidL;
        State UhatL = UL + gradX.col(elemL) * dxL(0) + gradY.col(elemL) * dxL(1);

        // Get state vector at the edge midpoint for right face
        const State& UR = U.col(elemR);
        Eigen::Vector2d centroidR = mesh.centroid(elemR);
        Eigen::Vector2d dxR = midpoint - centroidR;
        State UhatR = UR + gradX.col(elemR) * dxR(0) + gradY.col(elemR) * dxR(1);

        // use flux function to get the flux at the edge midpoint using uL and uR
        Eigen::Vector4d flux = computeNumericalFlux(UhatL, UhatR, normal);

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

        // REST NEEDS TO BE DONE
    }

    return R;
}

} // namespace solver