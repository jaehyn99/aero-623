// Calculates Barth-Jespersen limiter
#ifndef BJ_LIMITER_H
#define BJ_LIMITER_H

#include <Eigen/Dense>

class StateMesh;
std::vector<Eigen::Matrix<double,4,2>> BJLimiter(const std::vector<Eigen::Matrix<double,4,2>> Lgrad, const StateMesh& stateMesh);

#endif