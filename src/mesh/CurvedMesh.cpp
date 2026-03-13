#include "CurvedMesh.h"
#include "TriangularMesh.h"

#include "Eigen/Dense"

void CurvedMesh::curved_mesh(TriangularMesh& mesh, const Eigen::MatrixXi& B2E, int order)
{
    // loop over B2E and curve the edges and subsequently the nodes at the p=3 order
    for (int f_b = 0; f_b < B2E.rows(); f_b++)
    {
        // Compute the order of quadrature for internal points required as 2p+1 + 2*(Q-1)

        // Compute the order of quadrature for edge points required as 2p+1 + Q-1

        // get 1D Quadrature points along the reference edge (0,0)-(1,0) for the given order

        // get basis function values at quadrature points for each 1D lagrange polynomial

        // get the coordinates of the nodes from cubic spline on the edge and curve them using the basis functions

        // 

    }
}

void CurvedMesh::curved_mesh_alt(TriangularMesh& mesh, const Eigen::MatrixXi& B2E, int order)
{
    // loop over B2E and curve the edges at order p=3
    for (int f_b = 0; f_b < B2E.rows(); f_b++)
    {
        // get Jacobian of element mapping from reference element to physical element
        // get dx_vec/d(zeta) and dx_vec/d(eta) at the edge quadtrature points for base of reference traingle where d(zeta)/d(sig) = 1 and d(eta)/d(sig) = 0)

        // normalize this to get the unit tangent vector at the edge quadrature points

        // normal vector at the edge quadrature points is then given by rotating the tangent vector by 90 degrees
    }
}