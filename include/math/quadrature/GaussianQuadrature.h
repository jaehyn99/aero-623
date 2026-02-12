#ifndef GAUSSIAN_QUADRATURE_H
#define GAUSSIAN_QUADRATURE_H

#include "Quadrature.h"

template<int N = 1, int M = 1>
class GaussianQuadrature : public Quadrature<N, M>{
    public:
    using typename Quadrature<N, M>::DomainType;
    using typename Quadrature<N, M>::RangeType;
    using typename Quadrature<N, M>::Function;
    using typename Quadrature<N, M>::IntegrationDomain;

    GaussianQuadrature(std::size_t n): Quadrature<N, M>(), _n(n) {}
    virtual ~GaussianQuadrature() = default;

    virtual RangeType integrate(const Function& f, const IntegrationDomain& D) const override{
        if (!this->validIntegrationDomain(D)) throw std::invalid_argument("ERROR: Invalid integration domain.");

        RangeType sum;
        if constexpr (M == 1) sum = 0.0;
        else if constexpr (M == Eigen::Dynamic) sum = RangeType::Zero(M,1);
        else sum = RangeType::Zero();

        if constexpr (N == 1){
            auto nodes = getNodes(this->eval(D[0]), this->eval(D[1]));
            auto weights = getWeights(this->eval(D[0]), this->eval(D[1]));
            for (std::size_t i = 0; i < n(); i++) sum += weights[i] * f(nodes[i]);
        } else{
            constexpr std::size_t dim = (N != Eigen::Dynamic) ? N : D.size();
            std::size_t numNodes = 1; // total number of nodes in all dimensions;
            for (std::size_t i = 0; i < dim; i++) numNodes *= n();

            std::vector<DomainType> nodes(numNodes);
            if (N == Eigen::Dynamic){
                for (auto node: nodes) node.resize(dim);
            }
            std::vector<double> weights(numNodes, 1.0);
            this->getNodesAndWeights(0, dim, 0, numNodes, D, nodes, weights);
            for (std::size_t i = 0; i < numNodes; i++) sum += weights[i] * f(nodes[i]);
        }
        return sum;
    }

    protected:
    std::size_t _n; // number of nodes to evaluate
    virtual std::size_t n() const noexcept = 0; // number of actual nodes per dimension
    virtual Eigen::VectorXd getNodes(double l, double u) const = 0;
    virtual Eigen::VectorXd getWeights(double l, double u) const = 0;

    void getNodesAndWeights(std::size_t d, const std::size_t dim, std::size_t start, std::size_t stop,
                            const IntegrationDomain& D, std::vector<DomainType>& allNodes, std::vector<double>& allWeights) const{
        // Recursively form a Cartesian product matrix of nodes and their corresponding weights

        if (d == dim) return;

        double l = 0, u = 1;
        if (d == 0){
            l = this->eval(D[0]);
            u = this->eval(D[1]);
        } else if (d == 1){
            double param = allNodes[start][0];
            l = this->eval(D[2], param);
            u = this->eval(D[3], param);
        } else{
            Eigen::VectorXd param = allNodes[start](Eigen::seq(0, d));
            l = this->eval(D[2*d], param);
            u = this->eval(D[2*d+1], param);
        }
        if (std::isnan(l) || std::isnan(u)) throw std::runtime_error("ERROR: Cannot evaluate integration domain.");

        auto newNodes = getNodes(l, u);
        auto newWeights = getWeights(l, u);
        std::size_t n = std::round(std::pow(stop-start, 1.0/(dim-d)));
        std::size_t toFill = (stop-start)/n;
        for (auto i = start; i < stop; i++){
            // empty signifies either the nodes/weights list is not needed.
            allNodes[i][d] = newNodes[(i-start)/toFill];
            allWeights[i] *= newWeights[(i-start)/toFill];
        }

        for (std::size_t i = 0; i < n; i++)
            getNodesAndWeights(d+1, dim, start+i*toFill, start+(i+1)*toFill, D, allNodes, allWeights);
    }
};

#endif