#include "../include/ReadGRI.h"
#include "../include/ReadConnData.h"
#include "../include/MeshVerification.h"
#include <iostream>
#include "TriangularMesh.h"
#include <filesystem>

bool verificationSuite();

int main(){
    
    TriangularMesh mesh("projects/Project-1/mesh_coarse.gri");
    mesh.writeGri("projects/Project-1/mesh_coarse.gri");

    // std::cout << "C++ version: " << __cplusplus << std::endl;

    if (!verificationSuite()) throw std::runtime_error("Verification test(s) failed");
}

bool verificationSuite() {
    std::cout << "Running mesh verification suite.........................." << std::endl;
     // Mesh verification on test mesh

     std::filesystem::path currentDir = std::filesystem::current_path();

// Convert to a string and print
std::cout << "Current working directory: " << currentDir.string() << std::endl;

std::string testGriFile = "projects/Project-1/test.gri";
    std::vector<std::string> testTxtFiles = {"projects/Project-1/testperiodicEdges.txt",
            "projects/Project-1/testI2E.txt",
        "projects/Project-1/testIn.txt",
    "projects/Project-1/testB2E.txt",
"projects/Project-1/testBn.txt"};
    std::cout << "------------Mesh Verification for Test Grid------------" <<std::endl;
    meshVerification(testGriFile, testTxtFiles);
    std::cout <<std::endl;

// Mesh verification on coarse mesh
    std::string coarseGriFile = "projects/Project-1/mesh_coarse.gri";
    std::vector<std::string> coarseTxtFiles = {"projects/Project-1/mesh_coarseperiodicEdges.txt",
            "projects/Project-1/mesh_coarseI2E.txt",
        "projects/Project-1/mesh_coarseIn.txt",
    "projects/Project-1/mesh_coarseB2E.txt",
"projects/Project-1/mesh_coarseBn.txt"};
    std::cout << "------------Mesh Verification for Coarse Grid------------" <<std::endl;
    meshVerification(coarseGriFile, coarseTxtFiles);
    std::cout <<std::endl;

    // Mesh verification on refined mesh
    std::string localRefinedGriFile = "projects/Project-1/mesh_refined_2394.gri";
    std::vector<std::string> localRefinedTxtFiles = {"projects/Project-1/mesh_refined_2394periodicEdges.txt",
            "projects/Project-1/mesh_refined_2394I2E.txt",
        "projects/Project-1/mesh_refined_2394In.txt",
    "projects/Project-1/mesh_refined_2394B2E.txt",
"projects/Project-1/mesh_refined_2394Bn.txt"};
    std::cout << "------------Mesh Verification for Locally Refined Grid------------" <<std::endl;
    meshVerification(localRefinedGriFile, localRefinedTxtFiles);
    std::cout <<std::endl;

    std::string refined1GriFile = "projects/Project-1/meshGlobalRefined1.gri";
    std::vector<std::string> refined1TxtFiles = {"projects/Project-1/meshGlobalRefined1periodicEdges.txt",
            "projects/Project-1/meshGlobalRefined1I2E.txt",
        "projects/Project-1/meshGlobalRefined1In.txt",
    "projects/Project-1/meshGlobalRefined1B2E.txt",
"projects/Project-1/meshGlobalRefined1Bn.txt"};
    std::cout << "------------Mesh Verification for Global Refined Grid 1 ------------" <<std::endl;
    meshVerification(refined1GriFile, refined1TxtFiles);
    std::cout <<std::endl;

    std::string refined2GriFile = "projects/Project-1/meshGlobalRefined2.gri";
    std::vector<std::string> refined2TxtFiles = {"projects/Project-1/meshGlobalRefined2periodicEdges.txt",
            "projects/Project-1/meshGlobalRefined2I2E.txt",
        "projects/Project-1/meshGlobalRefined2In.txt",
    "projects/Project-1/meshGlobalRefined2B2E.txt",
"projects/Project-1/meshGlobalRefined2Bn.txt"};
    std::cout << "------------Mesh Verification for Global Refined Grid 2 ------------" <<std::endl;
    meshVerification(refined2GriFile, refined2TxtFiles);
    std::cout <<std::endl;

    std::string refined3GriFile = "projects/Project-1/meshGlobalRefined3.gri";
    std::vector<std::string> refined3TxtFiles = {"projects/Project-1/meshGlobalRefined3periodicEdges.txt",
            "projects/Project-1/meshGlobalRefined3I2E.txt",
        "projects/Project-1/meshGlobalRefined3In.txt",
    "projects/Project-1/meshGlobalRefined3B2E.txt",
"projects/Project-1/meshGlobalRefined3Bn.txt"};
    std::cout << "------------Mesh Verification for Global Refined Grid 3 ------------" <<std::endl;
    meshVerification(refined3GriFile, refined3TxtFiles);
    std::cout <<std::endl;

    return true;
}
