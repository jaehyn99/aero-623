#pragma once
#include <Eigen/Dense>

struct FluxResult {
    Eigen::Vector4d flux;
    double maxLambda;
};
