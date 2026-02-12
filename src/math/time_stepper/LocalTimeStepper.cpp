// #include "LocalTimeStepper.h"
// #include "StateMesh.h"
// #include "TriangularMesh.h"

// Eigen::ArrayXd LocalTimeStepper::dt(const StateMesh& u, const Eigen::ArrayXd& s) const noexcept{
//     auto mesh = u.mesh();
//     assert(mesh->numFaces() == s.size());
//     Eigen::ArrayXd steps = Eigen::ArrayXd::Zero(mesh->numElems());
//     const auto& edges = mesh->getFaces();
//     for (Eigen::Index i = 0; i < steps.size(); i++){
//         const TriangularMesh::Element& elem = mesh->elem(i);
//         double P = 0;
//         for (auto faceID: elem._faceID) P += s[faceID]*mesh->length(faceID);
//         double A = mesh->area(i);
//         steps[i] = 2*A*_minCFL/P;
//     }
//     return steps;
// }