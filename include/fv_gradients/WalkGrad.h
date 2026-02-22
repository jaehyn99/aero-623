// #pragma once

// #include "mesh/TriangularMesh.h"
// #include "mesh/StateMesh.h"
// #include "boundary_condition/BoundaryCondition.h"
// #include "FVGradient.h"


// #include <Eigen/Dense>
// #include <Eigen/StdVector>

// class WalkGrad : public FVGradient
// {
// public:

//     // Plane-normal gradient construction with hybrid treatment:
//     // - Interior and boundary faces: walk/Green-Gauss style accumulation using In

//     std::vector<Eigen::Matrix<double,4,2>> computeGradient(const Eigen::MatrixXi& I2E,
//                                                         const Eigen::MatrixXi& B2E,
//                                                         const Eigen::MatrixXd& In,
//                                                         const Eigen::MatrixXd& Bn,
//                                                         const Eigen::VectorXd& Area,
//                                                         const StateMesh& U) const override;
// };
