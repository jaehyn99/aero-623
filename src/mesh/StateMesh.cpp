#include "StateMesh.h"
#include "TriangularMesh.h"

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> mesh):
    _spatialMesh(mesh),
    _stateMesh(Eigen::ArrayXd::Zero(_spatialMesh->numElems()))
{}