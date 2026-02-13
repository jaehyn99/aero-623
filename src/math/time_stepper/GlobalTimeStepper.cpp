// #include "GlobalTimeStepper.h"
// #include "StateMesh.h"
// #include "TriangularMesh.h"

// Eigen::ArrayXd GlobalTimeStepper::dt(const StateMesh& u, const Eigen::ArrayXd& s) const noexcept{
//     auto mesh = u.mesh();
//     assert(mesh->numFaces() == s.size());
//     const auto& edges = mesh->getFaces();
//     double step = std::numeric_limits<double>::max();
//     for (auto elem: mesh->getElements()){
//         double P = 0;
//         for (auto faceID: elem._faceID) P += s[faceID]*mesh->length(faceID);
//         double A = elem._area;
//         step = std::min(step, 2*A*_minCFL/P);
//     }
//     return Eigen::ArrayXd::Constant(mesh->numElems(), step);
// }