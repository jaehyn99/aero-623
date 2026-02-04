// Reads connectivity data from .txt file
#pragma once
#include <vector>
#include <string>

class CONNData {
public:
    std::vector<std::vector<int>> I2E;
    std::vector<std::vector<double>> In;
    std::vector<double> sideLenInt;

    std::vector<std::vector<int>> B2E;
    std::vector<std::vector<double>> Bn;
    std::vector<double> sideLenBoundary;

    std::vector<int> periodicInd; //indices are 1-based
};