/**
 * Using adaptive Gauss-Kronrod Quadrature to evaluate the integral of f(x) over a finite domain.
**/

#ifndef GAUSS_KRONROD_H
#define GAUSS_KRONROD_H

#include "GaussLegendre.h"

class GaussKronrodBase : public GaussLegendreBase{
    // Dummy class to store the nodes and weights
    public:
    virtual ~GaussKronrodBase() = default;

    // Weights associated with the Gauss-Legendre nodes used in Gauss-Kronrod routines.
    inline static const std::vector<Eigen::ArrayXd> _Lweights =
        {{{0.8888888888888889}}, // starting at n = 1
         {{0.4909090909090909, 0.4909090909090909}},
         {{0.2684880898683338, 0.4509165386584746, 0.2684880898683338}},
         {{0.1700536053357231, 0.3269491896014521, 0.3269491896014521, 0.1700536053357231}},
         {{0.1152333166224727, 0.2410403392286477, 0.2829874178574913, 0.2410403392286477, 0.1152333166224727}}, // n = 5
         {{0.0836944404469058, 0.1810719943231379, 0.2337708641169940, 0.2337708641169940, 0.1810719943231379, 0.0836944404469058}},
         {{0.0630920926299783, 0.1406532597155260, 0.1903505780647856, 0.2094821410847274, 0.1903505780647856, 0.1406532597155260, 0.0630920926299783}},
         {{0.0494393950021392, 0.1116463708268402, 0.1566526061681888, 0.1814000250680342, 0.1814000250680342, 0.1566526061681888, 0.1116463708268402, 0.0494393950021392}},
         {{0.0396318951602617, 0.0907906816887264, 0.1300014068553415, 0.1564135277884832, 0.1648960128283498, 0.1564135277884832, 0.1300014068553415, 0.0907906816887264, 0.0396318951602617}},
         {{0.0325581623079644, 0.0750396748109197, 0.1093871588022977, 0.1347092173114738, 0.1477391049013377, 0.1477391049013377, 0.1347092173114738, 0.1093871588022977, 0.0750396748109197, 0.0325581623079644}}
        }; // n = 10 is the highest order

    // Extra Gauss-Kronrod nodes
    inline static const std::vector<Eigen::ArrayXd> _Knodes =
        {{{-0.7745966692414833,  0.7745966692414833}}, // starting at n = 1
         {{-0.9258200997725514,  0.0000000000000000,  0.9258200997725514}},
         {{-0.9604912687080204, -0.4342437493468024,  0.4342437493468024,  0.9604912687080204}},
         {{-0.9765602507375734, -0.6402862174963099,  0.0000000000000000,  0.6402862174963099,  0.9765602507375734}},
         {{-0.9840853600948420, -0.7541667265708494, -0.2796304131617834,  0.2796304131617834,  0.7541667265708494,  0.9840853600948420}}, // n = 5
         {{-0.9887032026126785, -0.8213733408650278, -0.4631182124753047,  0.0000000000000000,  0.4631182124753047,  0.8213733408650278,  0.9887032026126785}},
         {{-0.9914553711208126, -0.8648644233597690, -0.5860872354676913, -0.2077849550078988,  0.2077849550078988,  0.5860872354676913,  0.8648644233597690,  0.9914553711208126}},
         {{-0.9933798758817167, -0.8941209068474565, -0.6723540709451583, -0.3607010979281319,  0.0000000000000000,  0.3607010979281319,  0.6723540709451583,  0.8941209068474565,  0.9933798758817167}},
         {{-0.9946781606773405, -0.9149635072496780, -0.7344867651839337, -0.4754624791124598, -0.1642235636149869,  0.1642235636149869,  0.4754624791124598,  0.7344867651839337,  0.9149635072496780,  0.9946781606773405}},
         {{-0.9956571630258079, -0.9301574913557080, -0.7808177265864168, -0.5627571346686048, -0.2943928627014604,  0.0000000000000000,  0.2943928627014604,  0.5627571346686048,  0.7808177265864168,  0.9301574913557080,  0.9956571630258079}}
        }; // n = 10 is the highest order

    // Weights associated with the extra Gauss-Kronrod nodes
    inline static const std::vector<Eigen::ArrayXd> _Kweights =
        {{{0.5555555555555556, 0.5555555555555556}}, // starting at n = 1
         {{0.1979797979797980, 0.6222222222222222, 0.1979797979797980}},
         {{0.1046562260264674, 0.4013974147759618, 0.4013974147759618, 0.1046562260264674}},
         {{0.0629773736654729, 0.2667983404522846, 0.3464429818901361, 0.2667983404522846, 0.0629773736654729}}, 
         {{0.0425820367510821, 0.1868007965564924, 0.2728498019125588, 0.2728498019125588, 0.1868007965564924, 0.0425820367510821}}, // n = 5
         {{0.0303961541198201, 0.1373206046344471, 0.2132096522719623, 0.2410725801734653, 0.2132096522719623, 0.1373206046344471, 0.0303961541198201}},
         {{0.0229353220105295, 0.1047900103222504, 0.1690047266392678, 0.2044329400752993, 0.2044329400752993, 0.1690047266392678, 0.1047900103222504, 0.0229353220105295}},
         {{0.0178223833207102, 0.0824822989313582, 0.1362631092551723, 0.1720706085552109, 0.1844464057446914, 0.1720706085552109, 0.1362631092551723, 0.0824822989313582, 0.0178223833207102}},
         {{0.0143047756438387, 0.0665181559402740, 0.1117891346844181, 0.1452395883843660, 0.1628628274401153, 0.1628628274401153, 0.1452395883843660, 0.1117891346844181, 0.0665181559402740, 0.0143047756438387}},
         {{0.0116946388673721, 0.0547558965743521, 0.0931254545836974, 0.1234919762620657, 0.1427759385770600, 0.1494455540029175, 0.1427759385770600, 0.1234919762620657, 0.0931254545836974, 0.0547558965743521, 0.0116946388673721}}
        }; // n = 10 is the highest order
};

template<int N = 1, int M = 1>
class GaussKronrod : public GaussianQuadrature<N, M>, GaussKronrodBase{
    public:
    using typename Quadrature<N, M>::DomainType;
    using typename Quadrature<N, M>::RangeType;
    using typename Quadrature<N, M>::Function;
    using typename Quadrature<N, M>::IFunc1D;
    using typename Quadrature<N, M>::IFuncMD;
    using typename Quadrature<N, M>::IntegrationBound;
    using typename Quadrature<N, M>::IntegrationDomain;

    GaussKronrod(std::size_t n, double tol=1e-6): GaussianQuadrature<N, M>(n), GaussKronrodBase(), _tol(tol)
    {
        if (n < 1 || n > 10) throw std::domain_error("Error: Can only support 1 through 10 integration nodes.");
    }

    RangeType integrate(const Function& f, const IntegrationDomain& D) const override{
        if (!this->validIntegrationDomain(D)) throw std::invalid_argument("ERROR: Invalid integration domain.");
        return recursiveIntegral(f, D, _tol, 0);
    }

    protected:
    double _tol;

    std::size_t n() const noexcept override{ return 2*this->_n+1; }
    Eigen::VectorXd getNodes(double l, double u) const override{
        Eigen::ArrayXd x(n());
        x.head(this->_n) = _nodes[this->_n-1];
        x.tail(this->_n+1) = _Knodes[this->_n-1];
        return x * (u-l)/2 + (u+l)/2;
    }

    Eigen::VectorXd getWeights(double l, double u) const override{
        Eigen::ArrayXd w(n());
        w.head(this->_n) = _Lweights[this->_n-1];
        w.tail(this->_n+1) = _Kweights[this->_n-1];
        return w * (u-l)/2;
    }

    double error(const RangeType& GL, const RangeType& GK) const noexcept{
        // Empirical error estimation using Gauss-Legendre (GL)
        if constexpr (M == 1) return (GK == 0.0) ? std::abs(GL) : std::abs((GL-GK)/GK);
        else{
            RangeType diff = GL-GK;
            for (Eigen::Index i = 0; i < diff.size(); i++) diff[i] /= (GK[i] == 0) ? 1 : GK[i];
            return diff.norm();
        }
    }

    IntegrationBound midpoint(const IntegrationBound& l, const IntegrationBound& u) const{
        // Midpoint integration bound, used in adaptive steps
        if (std::get_if<double>(&l)){
            if (std::get_if<double>(&u)) return (this->eval(l) + this->eval(u)) / 2;
            if (std::get_if<IFunc1D>(&u)) return [&l, &u, this](double x) { return (this->eval(l) + this->eval(u,x)) / 2; };
            return [&l, &u, this](const Eigen::VectorXd& x) { return (this->eval(l) + this->eval(u,x)) / 2; };
        } else if (std::get_if<IFunc1D>(&l)) return [&l, &u, this](double x) { return (this->eval(l,x) + this->eval(u,x)) / 2; };
        return [&l, &u, this](const Eigen::VectorXd& x) { return (this->eval(l,x) + this->eval(u,x)) / 2; };
    }

    void recursiveSplitDomain(std::vector<IntegrationDomain>& newDomains, const IntegrationDomain& oldDomain,
                              const std::vector<IntegrationBound>& midpoints, std::size_t d, std::size_t start, std::size_t stop) const noexcept{
        if (d == midpoints.size()) return;
        auto pConst = std::get_if<double>(&midpoints[d]);
        if (pConst && std::isnan(*pConst)){
            // No split, midpoint stored as a nan
            for (std::size_t i = start; i < stop; i++){
                newDomains[i][2*d] = oldDomain[2*d];
                newDomains[i][2*d+1] = oldDomain[2*d+1];
            }
            recursiveSplitDomain(newDomains, oldDomain, midpoints, d+1, start, stop);
        } else{
            std::size_t diff = (stop-start)/2;
            for (std::size_t i = start; i < stop; i++){
                newDomains[i][2*d] = (i < start+diff) ? oldDomain[2*d] : midpoints[d];
                newDomains[i][2*d+1] = (i < start+diff) ? midpoints[d] : oldDomain[2*d+1];
            }
            recursiveSplitDomain(newDomains, oldDomain, midpoints, d+1, start, start+diff);
            recursiveSplitDomain(newDomains, oldDomain, midpoints, d+1, start+diff, stop);
        }
    }

    RangeType recursiveIntegral(const Function& f, const IntegrationDomain& D, double tol, std::size_t d) const{
        RangeType sum = GaussianQuadrature<N,M>::integrate(f, D);
        GaussLegendre<N,M> GL(this->_n);
        // would be nice to reuse repeated nodes, but in higher dimensions it doesn't help much
        double err = error(GL.integrate(f, D), sum);

        // Not yet converged, do adaptive quadrature
        if (err > tol){
            if constexpr (N == 1){
                double L = std::get<double>(D[0]);
                double U = std::get<double>(D[1]);
                std::size_t m = std::ceil(std::pow(err/tol, 1.0/(3*this->_n+1))); // number of subintervals
                double dx = (U-L)/m;
                sum *= 0; // reset the accumulator
                for (std::size_t i = 0; i < m; i++) sum += recursiveIntegral(f, {L+i*dx, L+(i+1)*dx}, tol/m, d);
            } else{
                constexpr std::size_t dim = (N != Eigen::Dynamic) ? N : D.size();
                std::size_t numSplits = std::ceil(std::log2(err/tol) / (3*this->_n+1)); // number of dimensions to split
                if (numSplits > dim) numSplits = dim;
                std::vector<IntegrationBound> midpoints(dim, std::numeric_limits<double>::quiet_NaN()); // stores the midpoints of the split dimensions
                for (std::size_t ns = 0; ns < numSplits; ns++){
                    midpoints[d] = midpoint(D[2*d], D[2*d+1]);
                    d = (d+1) % dim;
                }

                std::size_t m = std::pow(2, numSplits); // number of subdomains
                std::vector<IntegrationDomain> newDomains(m, IntegrationDomain(2*dim));
                recursiveSplitDomain(newDomains, D, midpoints, 0, 0, m);
                sum *= 0.0; // reset the accumulator
                for (auto domain: newDomains) sum += recursiveIntegral(f, domain, tol/m, d);
            }
        }
            
        return sum;
    };
};

#endif