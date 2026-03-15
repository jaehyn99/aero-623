#ifndef CURVED_MESH_H
#define CURVED_MESH_H

#include "Eigen/Dense"
#include <string>
#include <vector>
#include "TriangularMesh.h"

class CurvedMesh{
    public:
    void curved_mesh(TriangularMesh& mesh, const Eigen::MatrixXi& B2E, int order, int p);
};

#endif