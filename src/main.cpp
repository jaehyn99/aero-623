// src/main.cpp
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>

#include "TriangularMesh.h"

int main() {
    TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");
    mesh.writeGri("projects/Project-1/mesh_coarse.gri");

    return 0;
}