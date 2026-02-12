#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <Eigen/Dense>
#include <functional>
#include <variant>
#include <vector>

template<int N = 1, int M = 1,
        typename = std::enable_if_t<(N == Eigen::Dynamic || N >= 1) && (M == Eigen::Dynamic || M >= 1)>>
class Quadrature{
    public:
    using DomainType = std::conditional_t<N == 1, double, Eigen::Matrix<double, N, 1>>;
    using RangeType = std::conditional_t<M == 1, double, Eigen::Matrix<double, M, 1>>;
    using Function = std::function<RangeType(const DomainType&)>;

    using IFunc1D = std::function<RangeType(double)>; // 1D function as an integration bound
    using IFuncMD = std::function<RangeType(const Eigen::VectorXd&)>; // MD function as an integration bound
    using IntegrationBound = std::variant<RangeType, IFunc1D, IFuncMD>;
    using IntegrationDomain = std::vector<IntegrationBound>;

    virtual ~Quadrature() = default;
    virtual RangeType integrate(const Function&, const IntegrationDomain&) const = 0;

    protected:
    bool validIntegrationDomain(const IntegrationDomain& D) const noexcept{
        if (D.size() == 0 || D.size()%2 == 1) return false; // must contain a non-zero, even number of bounds
        if (N != Eigen::Dynamic && D.size() != N*2) return false; // the exact number of bounds is N*2 for fixed-size domains
        // First pair of bounds has to be constant
        if (!std::get_if<double>(&D[0]) || !std::get_if<double>(&D[1])) return false;
        // Second pair of bounds has to be constant or 1D function
        if (D.size() > 2){
            if (std::get_if<IFuncMD>(&D[2]) || std::get_if<IFuncMD>(&D[3])) return false;
        }
        // All bounds from then on have to be constant or MD function
        for (std::size_t i = 4; i < D.size(); i++){
            if (std::get_if<IFunc1D>(&D[i])) return false;            
        }
        return true;
    }

    RangeType eval(const IntegrationBound& B) const{
        if (auto pConst = std::get_if<RangeType>(&B)) return *pConst;
        throw std::bad_variant_access();
    }

    RangeType eval(const IntegrationBound& B, double param) const {
        if (auto pConst = std::get_if<RangeType>(&B)) return *pConst;
        if (auto pIFunc1D = std::get_if<IFunc1D>(&B)) return (*pIFunc1D)(param);
        throw std::bad_variant_access();
    }

    RangeType eval(const IntegrationBound& B, const Eigen::VectorXd& param) const {
        if (auto pConst = std::get_if<RangeType>(&B)) return *pConst;
        if (auto pIFuncMD = std::get_if<IFuncMD>(&B)) return (*pIFuncMD)(param);
        throw std::bad_variant_access();
    }
};

#endif