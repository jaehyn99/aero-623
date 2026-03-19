#include "StateMesh.h"
#include "TriangularMesh.h"
#include "FVGradient.h"

#include <fstream>

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, const Eigen::MatrixXd& array, Eigen::Index p): // std::shared_ptr<FVGradient> grad):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _p(p),
    _Np((p+1)*(p+2)/2),
    _stateMesh(array)
    //_gradientMethod(grad)
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols())/_Np);
    //computeGradient();
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::MatrixXd&& array, Eigen::Index p): //, std::shared_ptr<FVGradient> grad):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _p(p),
    _Np((p+1)*(p+2)/2),
    _stateMesh(std::move(array))
    //_gradientMethod(grad)
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols())/_Np);
    // computeGradient();
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::Index states, Eigen::Index p, double u0): //, std::shared_ptr<FVGradient> grad):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _p(p),
    _Np((p+1)*(p+2)/2),
    _stateMesh(Eigen::MatrixXd::Constant(states, _spatialMesh->numElems()*_Np, u0))
    //_gradientMethod(nullptr)
{} 

// std::vector<Eigen::Matrix<double,4,2>> StateMesh::computeGradient() const{
//     if (_gradientMethod) return _gradientMethod->computeGradient(*this);
//     throw std::invalid_argument("ERROR: No gradient computing function provided.");
// }


// std::vector<Eigen::Matrix<double,4,2>> StateMesh::computeGradient() const{
//     if (_gradientMethod)
//     {
//         // _spatialMesh->writeGri("placeholder_filename.gri");
//         // Eigen::MatrixXi I2E  = loadMatrix("placeholder_filenameI2E.txt", 4).cast<int>();
//         // Eigen::MatrixXi B2E  = loadMatrix("placeholder_filenameB2E.txt", 3).cast<int>();
//         // Eigen::MatrixXd In   = loadMatrix("placeholder_filenameIn.txt", 2);
//         // Eigen::MatrixXd Bn   = loadMatrix("placeholder_filenameBn.txt", 2);
//         // Eigen::VectorXd Area = loadVector("placeholder_filenameArea.txt");
//         return _gradientMethod->computeGradient(*this);
//     } throw std::invalid_argument("ERROR: No gradient computing function provided.");
// }