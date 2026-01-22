/**
 * Using Gauss-Legendre Quadrature to evaluate the integral of f(x) over a finite domain. 
**/

#ifndef GAUSS_LEGENDRE_H
#define GAUSS_LEGENDRE_H

#include "GaussianQuadrature.h"

#include <initializer_list>
#include <vector>

// Helper: construct a dynamic Eigen array from { ... }
inline Eigen::ArrayXd makeArrayXd(std::initializer_list<double> vals)
{
    Eigen::ArrayXd a(static_cast<Eigen::Index>(vals.size()));
    Eigen::Index i = 0;
    for (double v : vals) a(i++) = v;
    return a;
}

class GaussLegendreBase{
    // Dummy class to store the nodes and weights
    public:
    virtual ~GaussLegendreBase() = default;
    inline static const std::vector<Eigen::ArrayXd> _nodes =
        {
            makeArrayXd({ 0.0000000000000000 }), // starting at n = 1
            makeArrayXd({-0.5773502691896257,  0.5773502691896257}),
            makeArrayXd({-0.7745966692414833,  0.0000000000000000,  0.7745966692414833}),
            makeArrayXd({-0.8611363115940526, -0.3399810435848564,  0.3399810435848564,  0.8611363115940526}),
            makeArrayXd({-0.9061798459386639, -0.5384693101056832,  0.0000000000000000,  0.5384693101056832,  0.9061798459386639}), // n = 5
            makeArrayXd({-0.9324695142031517, -0.6612093864662645, -0.2386191860831968,  0.2386191860831968,  0.6612093864662645,  0.9324695142031517}),
            makeArrayXd({-0.9491079123427587, -0.7415311855993943, -0.4058451513773973,  0.0000000000000000,  0.4058451513773973,  0.7415311855993943,  0.9491079123427587}),
            makeArrayXd({-0.9602898564975364, -0.7966664774136267, -0.5255324099163293, -0.1834346424956497,  0.1834346424956497,  0.5255324099163293,  0.7966664774136267,  0.9602898564975364}),
            makeArrayXd({-0.9681602395076263, -0.8360311073266359, -0.6133714327005902, -0.3242534234038090,  0.0000000000000000,  0.3242534234038090,  0.6133714327005902,  0.8360311073266359,  0.9681602395076263}),
            makeArrayXd({-0.9739065285171715, -0.8650633666889844, -0.6794095682990245, -0.4333953941292472, -0.1488743389816312,  0.1488743389816312,  0.4333953941292472,  0.6794095682990245,  0.8650633666889844,  0.9739065285171715})
        }; // n = 10 is the highest order

    inline static const std::vector<Eigen::ArrayXd> _weights =
        {
            makeArrayXd({2.0000000000000000}), // starting at n = 1
            makeArrayXd({1.0000000000000000, 1.0000000000000000}),
            makeArrayXd({0.5555555555555556, 0.8888888888888889, 0.5555555555555556}),
            makeArrayXd({0.3478548451374543, 0.6521451548625463, 0.6521451548625463, 0.3478548451374543}),
            makeArrayXd({0.2369268850561888, 0.4786286704993663, 0.5688888888888889, 0.4786286704993663, 0.2369268850561888}), // n = 5
            makeArrayXd({0.1713244923791703, 0.3607615730481384, 0.4679139345726910, 0.4679139345726910, 0.3607615730481384, 0.1713244923791703}),
            makeArrayXd({0.1294849661688694, 0.2797053914892765, 0.3818300505051191, 0.4179591836734692, 0.3818300505051191, 0.2797053914892765, 0.1294849661688694}),
            makeArrayXd({0.1012285362903761, 0.2223810344533744, 0.3137066458778879, 0.3626837833783620, 0.3626837833783620, 0.3137066458778879, 0.2223810344533744, 0.1012285362903761}),
            makeArrayXd({0.0812743883615740, 0.1806481606948578, 0.2606106964029354, 0.3123470770400029, 0.3302393550012593, 0.3123470770400029, 0.2606106964029354, 0.1806481606948578, 0.0812743883615740}),
            makeArrayXd({0.0666713443086883, 0.1494513491505807, 0.2190863625159819, 0.2692667193099962, 0.2955242247147526, 0.2955242247147526, 0.2692667193099962, 0.2190863625159819, 0.1494513491505807, 0.0666713443086883})
        }; // n = 10 is the highest order

};

template<int N = 1, int M = 1>
class GaussLegendre : public GaussianQuadrature<N, M>, GaussLegendreBase{
    public:
    GaussLegendre(std::size_t n): GaussianQuadrature<N, M>(n), GaussLegendreBase()
    {
        if (n < 1 || n > 10) throw std::domain_error("Error: Can only support 1 through 10 integration nodes.");
    }

    protected:
    std::size_t n() const noexcept override{ return this->_n; }
    Eigen::VectorXd getNodes(double l, double u) const override{ return _nodes[this->_n-1] * (u-l)/2 + (u+l)/2; }
    Eigen::VectorXd getWeights(double l, double u) const override{ return _weights[this->_n-1] * (u-l)/2; }
};

#endif