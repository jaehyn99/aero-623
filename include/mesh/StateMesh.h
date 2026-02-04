// Array that stores the cell average of a state on each cell

#ifndef STATE_MESH_H
#define STATE_MESH_H

#include "Eigen/Dense"
#include <memory>

class TriangularMesh;
class StateMesh{
    public:
    StateMesh(std::shared_ptr<TriangularMesh> mesh); 

    double& operator[](Eigen::Index i) noexcept{ return _stateMesh[i]; }
    const double& operator[](Eigen::Index i) const noexcept{ return _stateMesh[i]; }

    protected:
    std::shared_ptr<TriangularMesh> _spatialMesh;
    Eigen::ArrayXd _stateMesh;
};

#endif