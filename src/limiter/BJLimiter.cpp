#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <memory>
// #include </mnt/c/Users/mmaru/Desktop/AE623/Project 2/external/eigen/Eigen/Dense>
// #include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/mesh/TriangularMesh.h"
// #include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/mesh/StateMesh.h"
// #include "/mnt/c/Users/mmaru/Desktop/AE623/Project 2/include/bj_limiter.h"
#include <Eigen/Dense>
#include "TriangularMesh.h"
#include "StateMesh.h"
#include "BJLimiter.h"


std::vector<Eigen::Matrix<double,4,2>> BJLimiter(const std::vector<Eigen::Matrix<double,4,2>> Lgrad, /*const TriangularMesh& triMesh,*/ const StateMesh& stateMesh) {
    // std::array<std::array<double>> Li(triMesh.numElems(),3);
    // For each cell:
    TriangularMesh& triMesh = *stateMesh.mesh();
    std::vector<Eigen::Matrix<double,4,2>> L_limit = Lgrad;
    for (std::size_t i=0;i<triMesh.numElems();i++){
    // Calculate u at each node
        // std::vector<std::vector<double>> uVals(4,4); //[u0, u1, u2, u3] at nodes; for each of the 4 states
        std::vector<std::vector<double>> uVals(4, std::vector<double>(4, 0)); //[u0, u1, u2, u3] at nodes; for each of the 4 states
        // std::vector<std::vector<double>> alphaVals(4,3);
        std::vector<std::vector<double>> alphaVals(4,std::vector<double>(3,0));
        Eigen::Vector2d centroid = triMesh.centroid(i);
        auto& elem = triMesh.elem(i);
        for (int k=0;k<4;k++){
            uVals[k][0] = stateMesh(k,i);//GET CURRENT STATE AT CELL i *******************
        }
        
        for (int j=0;j<3;j++) { // iterating over each of the nodes
        // Compute ray from centroid to node
            Eigen::Vector2d nodePoint = triMesh.node(elem._pointID[j]);
            Eigen::Vector2d rayVec = nodePoint - centroid;
            for (int k=0;k<4;k++){
                uVals[k][j+1] = uVals[k][0] + rayVec.dot(Lgrad[i].row(k).transpose()); // SOMEHOW GET L0 FOR THE CELL AND WANT TO TREAT LIKE AN EIGEN::Vector2d *******************
            }
        }

        std::vector<double> alphaSet;
        for (int k=0;k<4;k++) { // iterate over each state
            double umin, umax, alpha;
            // int idx;
            // idx = std::min_element(uVals[k].begin(),uVals[k].end());
            umin = *std::min_element(uVals[k].begin(),uVals[k].end()); // uVals[k][idx];
            // idx = std::max_element(uVals[k].begin(),uVals[k].end());
            umax = *std::max_element(uVals[k].begin(),uVals[k].end()); //uVals[k][idx];

            // find alpha
            for (int j=0;j<2;j++) { // iterate over 3 nodes
                if (uVals[k][j+1]-uVals[k][0] > 0) {
                    alphaVals[k][j] = std::min(1.0,(umax-uVals[k][0])/(uVals[k][j+1]-uVals[k][0]));
                }
                else if (uVals[k][j+1]-uVals[k][0] < 0) {
                    alphaVals[k][j] = std::min(1.0,(umin-uVals[k][0])/(uVals[k][j+1]-uVals[k][0]));
                }
                else {alphaVals[k][j] = 1;}
            }
            alpha = *std::min_element(alphaVals[k].begin(),alphaVals[k].end());
            L_limit[i](k) = Lgrad[i](k) * alpha;
        }
    }//repeat for all cells

    return L_limit;


};