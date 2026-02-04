#pragma once
#include <vector>
#include <string>

class GRIData {

public:
    struct Map {
        int nNode;
        int nElemTot;
        int Dim;
        int nBGroup;
        int nPG;
        std::vector<std::vector<double>> nodeXYZ; // Length: nNode
    };

    struct BoundaryGroup { // Length: nBGroup
        std::vector<int> nBFace;
        std::vector<int> nf;
        std::vector<std::string> Title;
        std::vector<std::vector<std::vector<int>>> NB;
    };

    struct ElementGroup { // Length: nElemGroup
        std::vector<int> nElem; // Number of elements in group i
        std::vector<int> order; // Geometry order of elements in egroup i (q)
        std::vector<std::string> basis; //String defining the geometry interpolation basis
        std::vector<std::vector<std::vector<int>>> NE; // List of nn(i) nodes for element j in egroup i
    };

    struct PeriodicGroup {
        std::vector<int> nPGNode;
        std::vector<std::string> periodicity;
        std::vector<std::vector<std::vector<int>>> NP; // Pairs for periodic nodes
    };

    struct Normals {
        std::vector<std::vector<double>> normalInterior;
        std::vector<std::vector<double>> normalExterior;
        std::vector<double> edgeLength;
    };

    Map map;
    BoundaryGroup boundaryGroup;
    ElementGroup elementGroup;
    PeriodicGroup periodicGroup;
    Normals normalGroup;

};