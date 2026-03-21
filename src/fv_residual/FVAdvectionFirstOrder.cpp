// #include "FVAdvectionFirstOrder.h"
// #include "InletBC.h"
// #include "InletOutletBC.h"
// #include "FVFlux.h"
// #include "StateMesh.h"
// #include "TriangularMesh.h"

// #include <iostream>
// #include <fstream>

// Eigen::MatrixXd FVAdvectionFirstOrder::computeResidual(const StateMesh& u) const{
//     auto mesh = u.mesh();
//     Eigen::MatrixXd residual = Eigen::MatrixXd::Zero(4, u.cellCount());
//     if (_flux){
//         std::vector<std::string> bcNames;
//         std::vector<bool> visited(mesh->numFaces(), false);

//         for (std::size_t i = 0; i < mesh->numFaces(); i++){
//             if (visited[i]) continue;
//             const auto& face = mesh->face(i);
//             std::size_t e = face._elemID[0];
//             Eigen::Vector2d normal = face._normal; // points outward, from elemL to elemR

//             if (face.isBoundaryFace() && !face.isPeriodicFace()){
//                 auto it = std::find(bcNames.cbegin(), bcNames.cend(), face._title);
//                 std::size_t boundaryID;
//                 if (it == bcNames.cend()){
//                     bcNames.push_back(face._title);
//                     boundaryID = bcNames.size() - 1;
//                 } else boundaryID = it - bcNames.cbegin();

//                 auto bc = u.bc(boundaryID);
//                 // Transient cases, these lines do nothing if running in steady-state
//                 if (auto inlet = dynamic_cast<InletBC*>(bc.get())){
//                     double y0 = mesh->node(face._pointID[0])[1];
//                     double y1 = mesh->node(face._pointID[1])[1];
//                     inlet->setTransientRho((y1+y0)/2);
//                 } else if (auto inlet = dynamic_cast<InletOutletBC*>(bc.get())){
//                     double y0 = mesh->node(face._pointID[0])[1];
//                     double y1 = mesh->node(face._pointID[1])[1];
//                     inlet->setTransientRho((y1+y0)/2);                       
//                 }
//                 residual.col(e) -= bc->computeFlux(u.cell(e), normal) * mesh->length(i) / mesh->area(e);
//             } else{
//                 std::size_t ne; // neighbor element index
//                 if (face.isPeriodicFace()){
//                     visited[i] = true;
//                     visited[face._periodicFaceID] = true;
//                     ne = face._periodicElemID;
//                 } else ne = face._elemID[1];

//                 Eigen::Vector4d flux = _flux->computeFlux(u.cell(e), u.cell(ne), normal);
//                 residual.col(e) -= flux * face._length / mesh->area(e);
//                 residual.col(ne) += flux * face._length / mesh->area(ne);
//             }
//         }
//     }
//     return residual;
// }