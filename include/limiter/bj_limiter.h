// Calculates Barth-Jespersen limiter
#pragma once
#include <vector>
#include <string>
#include <memory>

// #include </mnt/c/Users/mmaru/Desktop/AE623/Project 2/external/eigen/Eigen/Dense>
// #include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/mesh/TriangularMesh.h"
// #include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/mesh/StateMesh.h"
#include <Eigen/Dense>
#include "TriangularMesh.h"
#include "StateMesh.h"

std::vector<Eigen::Matrix<double,4,2>> bj_limiter(const std::vector<Eigen::Matrix<double,4,2>> Lgrad, /* const TriangularMesh& triMesh,*/ const StateMesh& stateMesh);