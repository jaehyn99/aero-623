#include "mesh/StateMesh.h"
#include "mesh/TriangularMesh.h"
#include <fstream>
#include "fv_gradients/Gradients.h"

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

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, const Eigen::MatrixXd& array):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _gradientMethod(nullptr),
    _stateMesh(array)
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols()));
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::MatrixXd&& array):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _gradientMethod(nullptr),
    _stateMesh(std::move(array))
{
    assert(_spatialMesh->numElems() == std::size_t(_stateMesh.cols()));
}

StateMesh::StateMesh(std::shared_ptr<TriangularMesh> spatialMesh, std::vector<std::shared_ptr<BoundaryCondition>>& bc, Eigen::Index states, double u0):
    _spatialMesh(spatialMesh),
    _bc(bc),
    _gradientMethod(nullptr),
    _stateMesh(Eigen::MatrixXd::Constant(states, _spatialMesh->numElems(), u0))
{} 

void StateMesh::computeGradients(){
    if (_gradientMethod)
    {
        _spatialMesh->writeGri("placeholder_filename.gri");
        Eigen::MatrixXi I2E  = loadMatrix("placeholder_filenameI2E.txt", 4).cast<int>();
        Eigen::MatrixXi B2E  = loadMatrix("placeholder_filenameB2E.txt", 3).cast<int>();
        Eigen::MatrixXd In   = loadMatrix("placeholder_filenameIn.txt", 2);
        Eigen::MatrixXd Bn   = loadMatrix("placeholder_filenameBn.txt", 2);
        Eigen::VectorXd Area = loadVector("placeholder_filenameArea.txt");
        _gradients = _gradientMethod->computeGradient(I2E, B2E, In, Bn, Area, *this);
    }
}