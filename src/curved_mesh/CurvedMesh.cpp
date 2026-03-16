// #include "CurvedMesh.h"
// #include "TriangularMesh.h"
// // #include "femQuadrature.hpp"
// // #include "Quadrature1D.hpp"
// #include "shape1D.h"
// #include "geom_snap.h"
// #include "Eigen/Dense"

// #include <iostream>

// void CurvedMesh::curved_mesh(const Eigen::MatrixXi& B2E, int Q, int p,
//                             const std::string& upperBladeFile,
//                             const std::string& lowerBladeFile,
//                             const std::string& upperCurveName,
//                             const std::string& lowerCurveName)
// {
//     // GET CURVE NODES FROM GEOM_SNAP ----------------------------------------------------
//     // --- Build splines from blade geometry files ---
//     auto upper = readBladePoints(upperBladeFile);
//     auto lower = readBladePoints(lowerBladeFile);
//     auto lowerShifted = offset(lower, 0.0, 18.0);

//     // Arc-length parameterization
//     auto sUpper = arcLength(upper);
//     auto sLower = arcLength(lowerShifted);

//     ParametricSpline2D upperSpline, lowerSpline;
//     upperSpline.sx.build(sUpper, getX(upper));
//     upperSpline.sy.build(sUpper, getY(upper));
//     lowerSpline.sx.build(sLower, getX(lowerShifted));
//     lowerSpline.sy.build(sLower, getY(lowerShifted));

//     // convert _nodes from Eigen::Vector2d to Point2D for geom_snap compatibility
//     std::vector<Point2D> meshNodes(_nodes.size());
//     for (std::size_t i = 0; i < _nodes.size(); i++) {
//         meshNodes[i] = {_nodes[i][0], _nodes[i][1]};
//     }

//     // build edge pairs from _faces directly instead of reading from gri
//     std::vector<BoundaryPair> upperPairs, lowerPairs;
//     for (const auto& face : _faces) {
//         if (face._title == upperCurveName)
//             upperPairs.push_back({face._pointID[0] + 1, face._pointID[1] + 1});
//         else if (face._title == lowerCurveName)
//             lowerPairs.push_back({face._pointID[0] + 1, face._pointID[1] + 1});
//     }

//     auto upperCurved = buildCurvedEdges(upperCurveName, upperPairs, meshNodes, 
//                                         upperSpline, sUpper.front(), sUpper.back(), Q);

//     auto lowerCurved = buildCurvedEdges(lowerCurveName, lowerPairs, meshNodes,
//                                         lowerSpline, sLower.front(), sLower.back(), Q);

    
//     // --- Combine into a single lookup: (n1,n2) -> curved nodes ---
//     // key: (n1, n2), value: vector of Q+1 Point2D
//     std::map<std::pair<int,int>, std::vector<Point2D>> curvedEdgeMap;
//     for (const auto& ce : upperCurved)
//         curvedEdgeMap[{ce.n1, ce.n2}] = ce.qnodes;
//     for (const auto& ce : lowerCurved)
//         curvedEdgeMap[{ce.n1, ce.n2}] = ce.qnodes;
//     // -----------------------------------------------------------------------------------


//     // CURVE THE MESH ---------------------------------------------------------------------

//     // Compute the order of quadrature for edge points required as 2p+1 + Q-1
//     int quad_order_edge = 2*p+1 + Q-1;

//     // Compute the order of quadrature for internal points required as 2p+1 + 2*(Q-1)
//     int quad_order_internal = 2*p+1 + 2*(Q-1);

//     // Initilize the 1D and 2D quadrature objects
//     Quadrature1D quad1D(quad_order_edge);
//     femQuadrature quad2D(quad_order_internal);

//     // 1D Quadrature points along the reference edge (0,0)-(1,0) for the given order
//     Eigen::VectorXd sigq_edge = quad1D.getQuadSig(quad_order_edge);
//     Eigen::VectorXd wq_edge = quad1D.getQuadW(quad_order_edge);
//     int nq_edge = wq_edge.size();

//     // 2D quadrature points and weights for the internal points for the given order
//     Eigen::MatrixXd xiq_internal = quad2D.getQuadXi(quad_order_internal);
//     Eigen::VectorXd wq_internal = quad2D.getQuadW(quad_order_internal);
//     int nq_internal = wq_internal.size();

//     // 1D basis functions and gradients of basis functions at quadrature points for each 1D lagrange polynomial
//     Eigen::MatrixXd phiq_edge = shapeL1D_quad(sigq_edge, Q);
//     Eigen::MatrixXd gphiq_edge = gradL1D_quad(sigq_edge, Q);

//     int nEdges = B2E.rows();
//     _curvedElems.reserve(nEdges);

//     // loop over B2E and curve the edges and subsequently the nodes at the Q=3 order
//     for (int f_b = 0; f_b < B2E.rows(); f_b++)
//     {
//         int bgroup = B2E(f_b, 2);
        
//         // only curve the blade edges, not the inlet and outlet
//         if (bgroup == 1 || bgroup == 5) 
//         {
//             // get element ID , local face ID, and global face ID for this boundary edge
//             int elemID = B2E(f_b, 0);
//             int localFaceID = B2E(f_b, 1);
//             int faceID = _elems[elemID]._faceID[localFaceID];

//             // get the two endpoints of this edge from the mesh nodes
//             int n1 = _faces[faceID]._pointID[0] + 1;
//             int n2 = _faces[faceID]._pointID[1] + 1;

//             // get physical curved node locations for this edge from the curvedEdgeMap
//             auto it = curvedEdgeMap.find({n1, n2});
//             if (it == curvedEdgeMap.end())
//                 it = curvedEdgeMap.find({n2, n1}); // try reversed order
            
//             if (it == curvedEdgeMap.end()) {
//                 std::cerr << "Error: No curved edge found for nodes (" << n1 << ", " << n2 << ")\n";
//                 continue;
//             }
//             const std::vector<Point2D>& qnodes = it->second;
//             Eigen::MatrixXd edge_nodes_curved(Q+1, 2);
//             for (int i = 0; i < Q+1; i++) {
//                 edge_nodes_curved(i, 0) = qnodes[i].x;
//                 edge_nodes_curved(i, 1) = qnodes[i].y;
//             }

//             // Physical coordinates at quadrature points
//             Eigen::MatrixXd xq_edge = phiq_edge * edge_nodes_curved;

//             // Tangent and normal vectors at edge quadrature points
//             Eigen::MatrixXd dx_dsig = gphiq_edge * edge_nodes_curved;
//             Eigen::MatrixXd tangents = dx_dsig.rowwise().normalized();
//             Eigen::MatrixXd normals = Eigen::MatrixXd(tangents.rows(), 2);
//             for (int i = 0; i < tangents.rows(); i++) 
//             {
//                 normals(i, 0) = -tangents(i, 1);
//                 normals(i, 1) = tangents(i, 0);
//             }

//             // Jacobian at edge quadrature points
//             Eigen::VectorXd edge_jacobians = dx_dsig.rowwise().norm();

//             // Jacobians at internal quadrature points
//             Eigen::MatrixXd xq_internal(nq_internal, 2);
            
//         }

//     }
// }