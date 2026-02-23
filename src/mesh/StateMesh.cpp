#include "StateMesh.h"
#include "TriangularMesh.h"
#include "FVGradient.h"

#include <fstream>

static Eigen::MatrixXd loadMatrix(const std::string& filename, int cols){
    std::ifstream file(filename);
    std::vector<double> values;
    double val;
    while (file >> val)
        values.push_back(val);
    int rows = values.size() / cols;
    Eigen::MatrixXd mat(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            mat(i,j) = values[i*cols + j];
    return mat;
}

static Eigen::VectorXd loadVector(const std::string& filename){
    std::ifstream file(filename);
    std::vector<double> values;
    double val;
    while (file >> val)
        values.push_back(val);
    Eigen::VectorXd vec(values.size());
    for (size_t i = 0; i < values.size(); ++i)
        vec(i) = values[i];
    return vec;
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, const Eigen::MatrixXd& array, std::shared_ptr<FVGradient> grad):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _stateMesh(array),
    _gradientMethod(grad)
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols()));
    computeGradient();
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::MatrixXd&& array, std::shared_ptr<FVGradient> grad):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _stateMesh(std::move(array)),
    _gradientMethod(grad)
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols()));
    computeGradient();
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::Index states, double u0, std::shared_ptr<FVGradient> grad):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _stateMesh(Eigen::MatrixXd::Constant(states, _spatialMesh->numElems(), u0)),
    _gradientMethod(nullptr)
{} 

// std::vector<Eigen::Matrix<double,4,2>> StateMesh::computeGradient() const{
//     if (_gradientMethod) return _gradientMethod->computeGradient(*this);
//     throw std::invalid_argument("ERROR: No gradient computing function provided.");
// }


std::vector<Eigen::Matrix<double,4,2>> StateMesh::computeGradient() const{
    if (_gradientMethod)
    {
        // _spatialMesh->writeGri("placeholder_filename.gri");
        // Eigen::MatrixXi I2E  = loadMatrix("placeholder_filenameI2E.txt", 4).cast<int>();
        // Eigen::MatrixXi B2E  = loadMatrix("placeholder_filenameB2E.txt", 3).cast<int>();
        // Eigen::MatrixXd In   = loadMatrix("placeholder_filenameIn.txt", 2);
        // Eigen::MatrixXd Bn   = loadMatrix("placeholder_filenameBn.txt", 2);
        // Eigen::VectorXd Area = loadVector("placeholder_filenameArea.txt");
        return _gradientMethod->computeGradient(*this);
    } throw std::invalid_argument("ERROR: No gradient computing function provided.");
}