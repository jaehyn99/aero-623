#ifndef STATE_MESH_H
#define STATE_MESH_H

#include "Eigen/Dense"

class BoundaryCondition;
class TriangularMesh;
class Gradients;

class StateMesh{
    public:
    // Constructs with the entire state matrix
    StateMesh(std::shared_ptr<TriangularMesh>, std::vector<std::shared_ptr<BoundaryCondition>>&, const Eigen::MatrixXd&);
    StateMesh(std::shared_ptr<TriangularMesh>, std::vector<std::shared_ptr<BoundaryCondition>>&, Eigen::MatrixXd&&);
    StateMesh(std::shared_ptr<TriangularMesh>, std::vector<std::shared_ptr<BoundaryCondition>>&, Eigen::Index=4, double=1.0);
    
    Eigen::Index size() const noexcept{ return _stateMesh.size(); }
    Eigen::Index cellCount() const noexcept{ return _stateMesh.cols(); }
    Eigen::Index stateCount() const noexcept{ return _stateMesh.rows(); }

    // Accessors and modifiers
    // Reference to a single value
    double& operator()(Eigen::Index s, Eigen::Index e) noexcept{ return _stateMesh(s,e); }
    const double& operator()(Eigen::Index s, Eigen::Index e) const noexcept{ return _stateMesh(s,e); }
    // View to a cell
    Eigen::MatrixXd::ColXpr cell(Eigen::Index e) noexcept{ return _stateMesh.col(e); }
    Eigen::MatrixXd::ConstColXpr cell(Eigen::Index e) const noexcept{ return _stateMesh.col(e); }
    // View to a state across all cells
    Eigen::MatrixXd::RowXpr state(Eigen::Index s) noexcept{ return _stateMesh.row(s); }
    Eigen::MatrixXd::ConstRowXpr state(Eigen::Index s) const noexcept{ return _stateMesh.row(s); }
    // Reference/View to the entire matrix
    Eigen::MatrixXd& matrix() noexcept{ return _stateMesh; }
    const Eigen::MatrixXd& matrix() const noexcept{ return _stateMesh; }
    auto flattened() noexcept{ return Eigen::Map<Eigen::VectorXd>(_stateMesh.data(), size()); }
    auto flattened() const noexcept{ return Eigen::Map<const Eigen ::VectorXd>(_stateMesh.data(), size()); }
    auto array() noexcept{ return _stateMesh.array(); }
    auto array() const noexcept{ return _stateMesh.array(); }
    // Reference to the underlying spatial mesh
    std::shared_ptr<TriangularMesh> mesh() const noexcept{ return _spatialMesh; }
    // Reference to the underlying boundary condition
    std::shared_ptr<BoundaryCondition> bc(std::size_t i) const noexcept{ return _bc[i]; }

    void setGradientMethod(std::shared_ptr<Gradients> gradientMethod){ _gradientMethod = gradientMethod; computeGradients(); }
    std::vector<Eigen::Matrix<double,4,2>>& getGradients() noexcept{ return _gradients; }
    const std::vector<Eigen::Matrix<double,4,2>>& getGradients() const noexcept{ return _gradients; }

    protected:
    std::shared_ptr<TriangularMesh> _spatialMesh;
    std::vector<std::shared_ptr<BoundaryCondition>> _bc; // should store the BCs in the order the non-periodic boundaries are listed in the input .GRI file
    Eigen::MatrixXd _stateMesh;
    std::vector<Eigen::Matrix<double,4,2>> _gradients; // precomputed gradients for the current state, to be used in second-order residuals
    std::shared_ptr<Gradients> _gradientMethod; // strategy for computing gradients, to be set by the user before calling computeGradients()
    void computeGradients(); // computes the gradients using the specified strategy and stores them in _gradients for later use in second-order residuals
};

#endif