#include "StateMesh.h"
#include "BoundaryCondition.h"
#include "TriangularMesh.h"

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, const Eigen::MatrixXd& array):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _stateMesh(array)
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols()));
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::MatrixXd&& array):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _stateMesh(std::move(array))
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols()));
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::Index states, double u0):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _stateMesh(Eigen::MatrixXd::Constant(states, _spatialMesh->numElems(), u0))
{}    