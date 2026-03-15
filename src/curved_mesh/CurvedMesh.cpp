#include "CurvedMesh.h"
#include "TriangularMesh.h"
#include "femQuadrature.hpp"
#include "Quadrature1D.hpp"
#include "shape1D.h"

#include "Eigen/Dense"

void CurvedMesh::curved_mesh(TriangularMesh& mesh, const Eigen::MatrixXi& B2E, int Q, int p)
{

    // Compute the order of quadrature for internal points required as 2p+1 + 2*(Q-1)
    int quad_order_internal = 2*p+1 + 2*(Q-1);

    // Compute the order of quadrature for edge points required as 2p+1 + Q-1
    int quad_order_edge = 2*p+1 + Q-1;

    // Initilize the 1D and 2D quadrature objects
    Quadrature1D quad1D(quad_order_edge);
    femQuadrature quad2D(quad_order_internal);

    // 1D Quadrature points along the reference edge (0,0)-(1,0) for the given order
    Eigen::VectorXd xiq_edge = quad1D.getQuadXi(quad_order_edge);
    Eigen::VectorXd wq_edge = quad1D.getQuadW(quad_order_edge);
    int nQ_edge = wq_edge.size();

    // 2D quadrature points and weights for the internal points for the given order
    Eigen::MatrixXd xiq_internal = quad2D.getQuadXi(quad_order_internal);
    Eigen::VectorXd wq_internal = quad2D.getQuadW(quad_order_internal);
    int nQ_internal = wq_internal.size();

    // 1D basis functions and gradients of basis functions at quadrature points for each 1D lagrange polynomial
    Eigen::MatrixXd phiq_edge = shapeL1D_quad(xiq_edge, Q);
    Eigen::MatrixXd gphiq_edge = gradL1D_quad(xiq_edge, Q);

    int nEdges = B2E.rows();
    mesh._curvedElems.reserve(nEdges);

    // loop over B2E and curve the edges and subsequently the nodes at the Q=3 order
    for (int f_b = 0; f_b < B2E.rows(); f_b++)
    {
        // get the coordinates of the nodes from cubic spline on the edge and curve them using the basis functions
        Eigen::MatrixXd edge_nodes_curved = 
        

        // 

    }
}

// void CurvedMesh::curved_mesh_alt(TriangularMesh& mesh, const Eigen::MatrixXi& B2E, int order)
// {
//     // loop over B2E and curve the edges at order p=3
//     for (int f_b = 0; f_b < B2E.rows(); f_b++)
//     {
//         // get Jacobian of element mapping from reference element to physical element
//         // get dx_vec/d(zeta) and dx_vec/d(eta) at the edge quadtrature points for base of reference traingle where d(zeta)/d(sig) = 1 and d(eta)/d(sig) = 0)

//         // normalize this to get the unit tangent vector at the edge quadrature points

//         // normal vector at the edge quadrature points is then given by rotating the tangent vector by 90 degrees
//     }
// }