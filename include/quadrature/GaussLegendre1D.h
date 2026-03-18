/**
 * Using Gauss-Legendre Quadrature to evaluate the integral of f(x) over a finite domain. 
**/

#ifndef GAUSS_LEGENDRE_1D_H
#define GAUSS_LEGENDRE_1D_H

#include "Quadrature.h"

template<int M = 1>
class GaussLegendre1D : public Quadrature<1, M>{
    public:
    using typename Quadrature<1, M>::RangeType;
    using typename Quadrature<1, M>::Function;
    using typename Quadrature<1, M>::IntegrationDomain;
    GaussLegendre1D(std::size_t ord): Quadrature<1, M>(), _ord(ord)
    {
        if (ord < 1 || ord > 13) throw std::domain_error("Error: Can only support up to a 13th-order of accuracy.");
    }

    RangeType integrate(const Function& f, const IntegrationDomain& D={0, 1}, const Eigen::VectorXd& J = Eigen::VectorXd(0)) const override{
        RangeType sum;
        if constexpr (M == 1) sum = 0.0;
        else if constexpr (M == Eigen::Dynamic) sum = RangeType::Zero(M,1);
        else sum = RangeType::Zero();

        auto nodes = getNodes(D);
        auto weights = getWeights(D);
        for (std::size_t i = 0; i < n(); i++){
            double Ji = J.size() == 0 ? 1 : J[i];
            sum += weights[i] * f(nodes[i]) * Ji;
        }
        return sum;
    }

    std::size_t n() const noexcept{ return std::ceil((_ord+1)/2.0); }
    std::size_t ord() const noexcept{ return _ord; }
    Eigen::VectorXd getNodes(const IntegrationDomain& D={0, 1}) const{ return _nodes[n()-1] * (D[1]-D[0]) + D[0]; }
    Eigen::VectorXd getWeights(const IntegrationDomain& D={0, 1}) const{ return _weights[n()-1] * (D[1]-D[0]); }

    protected:
    std::size_t _ord;
    inline static const std::vector<Eigen::ArrayXd> _nodes =
        {
            {{0.500000000000000}}, // starting at n = 1
            {{0.211324865405187, 0.788675134594813}},
            {{0.112701665379258, 0.500000000000000, 0.887298334620742}},
            {{0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026}},
            {{0.046910077030668, 0.230765344947158, 0.500000000000000, 0.769234655052841, 0.953089922969332}}, // n = 5
            {{0.033765242898424, 0.169395306766868, 0.380690406958402, 0.619309593041598, 0.830604693233132, 0.966234757101576}},
            {{0.025446043828621, 0.129234407200303, 0.297077424311301, 0.500000000000000, 0.702922575688699, 0.870765592799697, 0.974553956171379}}
        }; // n = 7 is the highest degrees, accurate up to order 13
    inline static const std::vector<Eigen::ArrayXd> _weights =
        {
            {{1.000000000000000}}, // starting at n = 1
            {{0.500000000000000, 0.500000000000000}},
            {{0.277777777777778, 0.444444444444444, 0.277777777777778}},
            {{0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727}},
            {{0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095}}, // n = 5
            {{0.085662246189585, 0.180380786524069, 0.233956967286345, 0.233956967286345, 0.180380786524069, 0.085662246189585}},
            {{0.064742483084435, 0.139852695744638, 0.190915025252560, 0.208979591836735, 0.190915025252560, 0.139852695744638, 0.064742483084435}}
        }; // n = 7 is the highest degrees, accurate up to order 13
};

#endif
