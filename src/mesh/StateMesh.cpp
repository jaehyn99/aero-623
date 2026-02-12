#include "StateMesh.h"
#include "TriangularMesh.h"

/* Storage order:
    ** Each column contains the states (mass, momentum, and energy) of the same element.
    ** Each row contains a single state across all elements.
*/

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> mesh):
    _spatialMesh(mesh),
    _stateMesh(Eigen::MatrixXd::Zero(4, _spatialMesh->numElems()))
{}

StateMesh& StateMesh::operator=(const Eigen::MatrixXd& other){
    assert(_stateMesh.rows() == other.rows() && _stateMesh.cols() == other.cols());
    _stateMesh = other;
    return *this;
}

StateMesh& StateMesh::operator=(Eigen::MatrixXd&& other){
    assert(_stateMesh.rows() == other.rows() && _stateMesh.cols() == other.cols());
    _stateMesh = std::move(other);
    return *this;
}