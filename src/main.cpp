#include "../include/ReadGRI.h"
#include "../include/ReadConnData.h"
#include "../include/project1Task3.h"
#include <iostream>

bool verificationSuite();

int main(){
    


    if (!verificationSuite()) throw std::runtime_error("Verification test(s) failed");
}

bool verificationSuite() {
    std::cout << "Running mesh verification suite.........................." << std::endl;
     // Mesh verification on test mesh
    std::string testGriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/test.gri";
    std::vector<std::string> testTxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testperiodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testI2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testIn.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testB2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/testBn.txt"};
    std::cout << "------------Mesh Verification for Test Grid------------" <<std::endl;
    meshVerification(testGriFile, testTxtFiles);
    std::cout <<std::endl;

// Mesh verification on coarse mesh
    std::string coarseGriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_coarse.gri";
    std::vector<std::string> coarseTxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_coarseperiodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_coarseI2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_coarseIn.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_coarseB2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_coarseBn.txt"};
    std::cout << "------------Mesh Verification for Coarse Grid------------" <<std::endl;
    meshVerification(coarseGriFile, coarseTxtFiles);
    std::cout <<std::endl;

    // Mesh verification on refined mesh
    std::string localRefinedGriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394.gri";
    std::vector<std::string> localRefinedTxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394periodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394I2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394In.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394B2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/mesh_refined_2394Bn.txt"};
    std::cout << "------------Mesh Verification for Locally Refined Grid------------" <<std::endl;
    meshVerification(localRefinedGriFile, localRefinedTxtFiles);
    std::cout <<std::endl;

    std::string refined1GriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined1.gri";
    std::vector<std::string> refined1TxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined1periodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined1I2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined1In.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined1B2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined1Bn.txt"};
    std::cout << "------------Mesh Verification for Global Refined Grid 1 ------------" <<std::endl;
    meshVerification(refined1GriFile, refined1TxtFiles);
    std::cout <<std::endl;

    std::string refined2GriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined2.gri";
    std::vector<std::string> refined2TxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined2periodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined2I2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined2In.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined2B2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined2Bn.txt"};
    std::cout << "------------Mesh Verification for Global Refined Grid 2 ------------" <<std::endl;
    meshVerification(refined2GriFile, refined2TxtFiles);
    std::cout <<std::endl;

    std::string refined3GriFile = "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined3.gri";
    std::vector<std::string> refined3TxtFiles = {"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined3periodicEdges.txt",
            "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined3I2E.txt",
        "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined3In.txt",
    "/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined3B2E.txt",
"/mnt/c/Users/mmaru/Desktop/AE623/Project 1/projects/Project-1/meshGlobalRefined3Bn.txt"};
    std::cout << "------------Mesh Verification for Global Refined Grid 3 ------------" <<std::endl;
    meshVerification(refined3GriFile, refined3TxtFiles);
    std::cout <<std::endl;

    return true;
}