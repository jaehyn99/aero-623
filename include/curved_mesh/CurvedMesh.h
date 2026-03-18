#ifndef CURVED_MESH_H
#define CURVED_MESH_H

#include "TriangularMesh.h"
#include "Eigen/Dense"
#include <vector>
#include <iostream>

class CurvedMesh : public TriangularMesh {
public:
    CurvedMesh(const std::string& file_name) : TriangularMesh(file_name) {}
    
    void CurvedMesh::curved_mesh(const Eigen::MatrixXi& B2E, int Q, int p,
                                const std::string& upperBladeFile,
                                const std::string& lowerBladeFile,
                                const std::string& upperCurveName,
                                const std::string& lowerCurveName);
    };

#endif