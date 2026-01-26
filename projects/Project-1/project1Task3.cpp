// #include "../../include/ReadGRI.h"
#include "../include/ReadGRI.h"
#include <fstream>
#include <sstream>
#include <iostream>

GRIData readGriFile(const std::string& filename) {

    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Could not open file " + filename);

    Map gridMap;

    file >> gridMap.nNode >> gridMap.nElemTot >> gridMap.Dim;

    // Saves node coordinates
    std::vector<double> point;
    for (int i=0;i<gridMap.nNode;i++){
        double x, y, z = 0.0;
        file >> x >> y;
        if (gridMap.Dim == 3) {
            file >> z;
        }
        point = {x, y, z};
        gridMap.nodeXYZ.push_back(point);
    }

    // Saving boundary data
    file >> gridMap.nBGroup;
    int nBFace, nf;
    std::string title;
    std::vector<int> nb;
    int nbTemp;
    BoundaryGroup boundary;
    for (int i=0;i<gridMap.nBGroup;i++) {
        file >> nBFace >> nf >> title;
        // std::string line;
        // if (std::getline(file, line)) {  // Read one line from file
        //     std::istringstream iss(line); // Use stringstream to parse the line
        //     iss >> nBFace >> nf >> title; // Extract the first two ints and the string
        // }
        for (int j=0;j<nBFace;j++) {
            file >> nbTemp;
            nb.push_back(nbTemp);
        }
        boundary


    }
}

// ReadGRI::readFile(const std::string& filename)
//     : filename_(filename) {}

// bool ReadGRI::readFile() {
//     std::ifstream file(filename_);
//     if (!file.is_open()) {
//         std::cerr << "Error opening file: " << filename_ << std::endl;
//         return false;
//     }

