#include "../include/ReadGRI.h"
#include "../include/ReadConnData.h"
#include "../include/MeshVerification.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

GRIData readGriFile(const std::string& filename) {

    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Could not open file " + filename);

    GRIData data;

    file >> data.map.nNode >> data.map.nElemTot >> data.map.Dim;

    // Saves node coordinates
    data.map.nodeXYZ.reserve(data.map.nNode);
    for (int i=0;i<data.map.nNode;i++){
        double x, y, z = 0.0;
        file >> x >> y;
        if (data.map.Dim == 3) {
            file >> z;
        }
        data.map.nodeXYZ.emplace_back(std::vector<double>{x, y, z});
    }

    // Saving boundary data
    file >> data.map.nBGroup;
    data.boundaryGroup.nBFace.reserve(data.map.nBGroup);
    data.boundaryGroup.nf.reserve(data.map.nBGroup);
    data.boundaryGroup.Title.reserve(data.map.nBGroup);
    data.boundaryGroup.NB.reserve(data.map.nBGroup);
    for (int i=0;i<data.map.nBGroup;i++) {
        int nBFaceTemp, nfTemp;
        std::string title;
        file >> nBFaceTemp >> nfTemp >> title;

        std::vector<std::vector<int>> nbVec;
        nbVec.reserve(nBFaceTemp);
        for (int j=0;j<nBFaceTemp;j++) {
            std::vector<int> innerVec;
            innerVec.reserve(nfTemp);
            for (int k=0;k<nfTemp;k++) {
                int nbTemp;
                file >> nbTemp;
                innerVec.emplace_back(nbTemp);
            }
            nbVec.emplace_back(innerVec);
        }
        data.boundaryGroup.NB.emplace_back(std::move(nbVec));
        data.boundaryGroup.nBFace.emplace_back(nBFaceTemp);
        data.boundaryGroup.nf.emplace_back(nfTemp);
        data.boundaryGroup.Title.emplace_back(std::move(title));
    }

    // Saving element data
    int nElemCount = 0;
    while (nElemCount < data.map.nElemTot) {
        int nElemTemp, orderTemp;
        std::string basisTemp;
        file >> nElemTemp >> orderTemp >> basisTemp;
        std::vector<std::vector<int>> neVec;
        neVec.reserve(nElemTemp);
        nElemCount = nElemCount+nElemTemp;
        for (int j=0;j<nElemTemp;j++) {
            std::vector<int> innerVec2;
            int nn = 3; // CAUTION: This is hard-coded for triangular mesh right now!!!!!
            innerVec2.reserve(nn);
            for (int k=0;k<nn;k++) {
                int neTemp;
                file >> neTemp;
                innerVec2.emplace_back(neTemp);
            }
            neVec.emplace_back(innerVec2);
        }
        data.elementGroup.nElem.emplace_back(nElemTemp);
        data.elementGroup.order.emplace_back(orderTemp);
        data.elementGroup.basis.emplace_back(std::move(basisTemp));
        data.elementGroup.NE.emplace_back(std::move(neVec));
    };

    // Saving periodicity data
    
    std::string declarePeriodicity;
    file >> data.map.nPG >> declarePeriodicity;
    if (declarePeriodicity == "PeriodicGroup") {
        data.periodicGroup.nPGNode.reserve(data.map.nPG);
        data.periodicGroup.periodicity.reserve(data.map.nPG);
        data.periodicGroup.NP.reserve(data.map.nPG);
        for (int i=0;i<data.map.nPG;i++) {
            int npgNodeTemp;
            std::string periodTemp;
            std::vector<std::vector<int>> innerVec3;
            file >> npgNodeTemp >> periodTemp;
            for (int j=0;j<npgNodeTemp;j++) {
                int pair1, pair2;
                file >> pair1 >> pair2;
                innerVec3.emplace_back(std::vector<int>{pair1,pair2});
            }
            data.periodicGroup.NP.emplace_back(std::move(innerVec3));
            data.periodicGroup.nPGNode.emplace_back(npgNodeTemp);
            data.periodicGroup.periodicity.emplace_back(std::move(periodTemp));
        }
    }
    else {data.map.nPG = -1;}
    // std::cout << "saved GRI data successfully " << std::endl;
    return data;

}



CONNData readConnData(const std::vector<std::string>& filenames) {
// File names are stored as: {"testperiodicEdges.txt", "I2E.txt", "In.txt", "B2E.txt", "Bn.txt"}
    CONNData data;
    // Saving periodic indices
    std::ifstream file(filenames[0]);
    if (!file) throw std::runtime_error("Could not open file " + filenames[0]);
    std::string line;
    int count = 0;
    while(std::getline(file,line)) {
        count ++;
    }
    file.clear();
    file.seekg(0);
    data.periodicInd.resize(count);
    for (int i=0;i<count;i++) {
        file >> data.periodicInd[i];
    }
    file.close();
    
    // Saving I2E data
    std::ifstream file1(filenames[1]);
    if (!file1) throw std::runtime_error("Could not open file " + filenames[1]);
    int nElemInt = 0;
    while(std::getline(file1,line)) {
        nElemInt ++;
    }
    file1.clear();
    file1.seekg(0);
    int realElemInt = 0;
    data.I2E.reserve(nElemInt);
    for (int i=0;i<nElemInt;i++) {
        // Read in I2E file
        std::vector<int> tempI2E(4);
        for (int j=0;j<4;j++) {
            file1 >> tempI2E[j];
        }
        if (i<count) {
            if (i<data.periodicInd[i]) { // periodicInd is 1-based
                data.I2E.emplace_back(tempI2E);
                realElemInt ++;
            }
        }
        else {
            data.I2E.emplace_back(tempI2E);
            realElemInt ++;
        }
    }
    file1.close();


    // Saving In data
    std::ifstream file2(filenames[2]);
    if (!file2) throw std::runtime_error("Could not open file " + filenames[2]);
    data.In.reserve(nElemInt);
    for (int i=0;i<nElemInt;i++) {
        // Read in I2E file
        double point1, point2;
        file2 >> point1 >> point2;
        if (i<count) {
            if (i<data.periodicInd[i]) { // periodicInd is 1-based
                data.In.emplace_back(std::vector<double> {point1, point2});
            }
        }
        else {
            data.In.emplace_back(std::vector<double> {point1, point2});
        }
    }
    file2.close();

    // Saving B2E data
    std::ifstream file3(filenames[3]);
    if (!file3) throw std::runtime_error("Could not open file " + filenames[3]);
    int nElemBoundary = 0;
    while(std::getline(file3,line)) {
        nElemBoundary ++;
    }
    file3.clear();
    file3.seekg(0);
    data.B2E.reserve(nElemBoundary);
    for (int i=0;i<nElemBoundary;i++) {
        // Read in I2E file
        std::vector<int> tempB2E(4);
        for (int j=0;j<3;j++) {
            file3 >> tempB2E[j];
        }
        data.B2E.emplace_back(tempB2E);
    }
    file3.close();

    // Saving Bn data
    std::ifstream file4(filenames[4]);
    if (!file4) throw std::runtime_error("Could not open file " + filenames[4]);
    data.Bn.reserve(nElemBoundary);
    for (int i=0;i<nElemBoundary;i++) {
        // Read in I2E file
        double point1, point2;
        file4 >> point1 >> point2;
        data.Bn.emplace_back(std::vector<double> {point1, point2});
    }
    file4.close();
    // std::cout << "Successfully read conn data" << std::endl;
    return data;
}

bool meshVerification(const std::string& GriFile, const std::vector<std::string>& txtFiles) {
    // std::cout << "Running mesh verification"<<std::endl;
    
    GRIData gridData = readGriFile(GriFile);
    CONNData connData = readConnData(txtFiles);
    std::cout << "Number of elements: " << gridData.map.nElemTot << std::endl;

    // Calculate lengths of sides for boundary and interior edges
    int lenInt = connData.In.size();
    int lenBoun = connData.Bn.size();
    connData.sideLenInt.reserve(lenInt);
    // std::cout << "Reserved side length array: "<< lenInt << std::endl;
    for (int i=0;i<lenInt;i++) {
        // std::cout << "iter: " << i << std::endl;
        int elemL = connData.I2E[i][0];
        int faceL = connData.I2E[i][1];
        // std::cout << "determined elem face" << std::endl;
        // Note: Nodes are 1-based
        // int indAcross = gridData.elementGroup.NE[0][elemL-1][faceL-1];
        int ind1, ind2;
        if (faceL == 1) {
            ind1 = gridData.elementGroup.NE[0][elemL-1][1];
            ind2 = gridData.elementGroup.NE[0][elemL-1][2];
        }
        else if (faceL == 2 ){
            ind1 = gridData.elementGroup.NE[0][elemL-1][0];
            ind2 = gridData.elementGroup.NE[0][elemL-1][2];
        }
        else if (faceL == 3) {
            ind1 = gridData.elementGroup.NE[0][elemL-1][0];
            ind2 = gridData.elementGroup.NE[0][elemL-1][1];
        }
        else {throw std::runtime_error("Face index is not stored as 1-2-3 ");};
        // std::cout << "determined ind1 ind2: "<< ind1<< " " << ind2 << std::endl;
        // std::cout << gridData.map.nodeXYZ[ind1-1][0] << std::endl;
        // std::cout << gridData.map.nodeXYZ[ind1-1][1] << std::endl;
        // std::cout << gridData.map.nodeXYZ[ind2-1][0] << std::endl;
        // std::cout << gridData.map.nodeXYZ[ind2-1][1] << std::endl;
        double lenSeg = sqrt(pow(gridData.map.nodeXYZ[ind2-1][0]-gridData.map.nodeXYZ[ind1-1][0],2) + pow(gridData.map.nodeXYZ[ind2-1][1]-gridData.map.nodeXYZ[ind1-1][1],2));
        // std::cout << "calculated length" << std::endl;
        connData.sideLenInt.emplace_back(lenSeg);
        // std::cout << "saved side length" << std::endl;
    }
    // std::cout << "successfully computed lengths of interior sides" << std::endl;

    connData.sideLenBoundary.reserve(lenBoun);
    for (int i=0;i<lenBoun;i++) {
        int elem = connData.B2E[i][0]; // element adjacent to boundary face
        int face = connData.B2E[i][1]; // boundary face's local face number within that element
        // int bgroup = connData.B2E[i][2]; // boundary index of the face

        // Note: Nodes are 1-based
        // int indAcross = gridData.elementGroup.NE[0][elem-1][face-1];
        int ind1, ind2;

        if (face == 1) {
            ind1 = gridData.elementGroup.NE[0][elem-1][1];
            ind2 = gridData.elementGroup.NE[0][elem-1][2];
        }
        else if (face == 2 ){
            ind1 = gridData.elementGroup.NE[0][elem-1][0];
            ind2 = gridData.elementGroup.NE[0][elem-1][2];
        }
        else if (face == 3) {
            ind1 = gridData.elementGroup.NE[0][elem-1][0];
            ind2 = gridData.elementGroup.NE[0][elem-1][1];
        }
        else {std::cout << "faceL = " << face << std::endl; throw std::runtime_error("Face index is not stored as 1-2-3 ");};

        double lenSeg = sqrt(pow(gridData.map.nodeXYZ[ind2-1][0]-gridData.map.nodeXYZ[ind1-1][0],2) + pow(gridData.map.nodeXYZ[ind2-1][1]-gridData.map.nodeXYZ[ind1-1][1],2));
        connData.sideLenBoundary.emplace_back(lenSeg);
    }
    // std::cout << "successfully computed lengths of boundary sides" << std::endl;

    // Calculate error in each grid cell
    std::vector<double> Ee_vecX(gridData.map.nElemTot,0.0);
    std::vector<double> Ee_vecY(gridData.map.nElemTot,0.0);
    std::vector<double> Ee;
    Ee.reserve(gridData.map.nElemTot);
    // Loop over interior edges (boarder two cells each)
    for (int i=0;i<lenInt;i++) {
        double errLeftX, errLeftY, errRightX, errRightY;
        errRightX = -1*connData.sideLenInt[i] * connData.In[i][0];
        errRightY = -1*connData.sideLenInt[i] * connData.In[i][1];
        errLeftX = connData.sideLenInt[i] * connData.In[i][0];
        errLeftY = connData.sideLenInt[i] * connData.In[i][1];
        Ee_vecX[connData.I2E[i][2]-1] = Ee_vecX[connData.I2E[i][2]-1]+errRightX;
        Ee_vecY[connData.I2E[i][2]-1] = Ee_vecY[connData.I2E[i][2]-1]+errRightY;
        Ee_vecX[connData.I2E[i][0]-1] = Ee_vecX[connData.I2E[i][0]-1]+errLeftX;
        Ee_vecY[connData.I2E[i][0]-1] = Ee_vecY[connData.I2E[i][0]-1]+errLeftY;
    }
    // Loop over exterior edges (only one cell)
    for (int i=0;i<lenBoun;i++) {
        // std::cout << "iter: " << i << std::endl; // seg fault at 296
        double errX, errY;
        errX = connData.sideLenBoundary[i] * connData.Bn[i][0];
        errY = connData.sideLenBoundary[i] * connData.Bn[i][1];
        Ee_vecX[connData.B2E[i][0]-1] = Ee_vecX[connData.B2E[i][0]-1]+errX;
        Ee_vecY[connData.B2E[i][0]-1] = Ee_vecY[connData.B2E[i][0]-1]+errY;
    }
    // std::cout << "successfully computed cell errors" << std::endl;

    double maxErr = 0.0;
    int indMaxErr = -1;
    for (int i=0;i<gridData.map.nElemTot;i++) {
        Ee[i] = sqrt(pow(Ee_vecX[i],2) + pow(Ee_vecY[i],2));
        // std::cout << "Ex = " << Ee_vecX[i] << " Ey = " << Ee_vecY[i] << " Ee = " << Ee[i] << std::endl;
        if (Ee[i] > maxErr) {
            maxErr = Ee[i];
            indMaxErr = i+1; // 1-based indexing for cell numbering
        }
    }
    std::cout << "Maximum error in cell #: " << indMaxErr << "  Error: " << Ee[indMaxErr-1] << std::endl;
    return true;

}


