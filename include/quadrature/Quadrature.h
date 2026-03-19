#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <Eigen/Dense>
#include <functional>
#include <vector>

template<int N = 1, int M = 1,
        typename = std::enable_if_t<(N == 1 || N == 2) && (M == Eigen::Dynamic || M >= 1)>>
class Quadrature{
    public:
    using RangeType = std::conditional_t<M == 1, double, Eigen::Matrix<double, M, 1>>;
    using Function = std::conditional_t<N == 1, std::function<RangeType(double)>, std::function<RangeType(double, double)>>;
    using IntegrationDomain = Eigen::Matrix<double, 2*N, 1>; // Vector2d if N = 1, Vector4d if N = 2

    virtual ~Quadrature() = default;
    // J is the vector of Jacobian determinants at the corresponding quadrature nodes, used to scale the contribution of each quadrature node
    virtual RangeType integrate(const Function&, const IntegrationDomain&, const Eigen::VectorXd& = Eigen::VectorXd(0)) const = 0;
};

#endif