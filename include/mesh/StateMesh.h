// Array that stores the cell averages of a state on each cell

#ifndef STATE_MESH_H
#define STATE_MESH_H

#include "Eigen/Dense"
#include <memory>

class TriangularMesh;
class StateMesh{
    public:
    StateMesh(std::shared_ptr<TriangularMesh> mesh); 

    double size() const noexcept{ return _stateMesh.size(); }
    double& operator()(Eigen::Index var, Eigen::Index elem) noexcept{ return _stateMesh(var, elem); }
    const double& operator()(Eigen::Index var, Eigen::Index elem) const noexcept{ return _stateMesh(var, elem); }
    Eigen::MatrixXd& matrix() noexcept{ return _stateMesh; }
    const Eigen::MatrixXd& matrix() const noexcept{ return _stateMesh; }

    StateMesh& operator=(const Eigen::MatrixXd& other);
    StateMesh& operator=(Eigen::MatrixXd&& other);
    std::shared_ptr<TriangularMesh> mesh() const noexcept{ return _spatialMesh; }

    protected:
    std::shared_ptr<TriangularMesh> _spatialMesh;
    Eigen::MatrixXd _stateMesh;
};

#endif