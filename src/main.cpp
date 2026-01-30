#include "TriangularMesh.h"
#include <iostream>

int main(){
    TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");
    mesh.writeGri("projects/Project-1/mesh_coarse.gri");
    
    // std::cout << "Node points:" << std::endl;
    // for (auto node: mesh.getNodes()) std::cout << node.transpose() << std::endl;

    // std::cout << std::endl;
    // std::cout << "Faces" << std::endl;
    // for (auto face: mesh.getFaces()){
    //     std::cout << "Face containing points " << face._pointID[0] << " and " << face._pointID[1];
    //     if (face.isBoundaryFace()){
    //         if (face._periodicFaceID == -1) std::cout << ", and is a boundary face";
    //         else std::cout << ", and is a periodic boundary face";
    //     }
    //     std::cout << ". It borders element(s) " << face._elemID[0] << " and " << face._elemID[1] << "." << std::endl;
    // }

    // std::cout << std::endl;
    // std::cout << "Elements" << std::endl;
    // for (auto elem: mesh.getElements()){
    //     std::cout << "Element containing points " << elem._pointID[0] << ", " << elem._pointID[1] << " and " << elem._pointID[2] << ". ";
    //     std::cout << "It contains faces " << elem._faceID[0] << ", " << elem._faceID[1] << " and " << elem._faceID[2] << "." << std::endl;
    // }
}