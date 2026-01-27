#include "mesh/GriIO.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace mesh {

GriMesh2D readGri2D(const std::string& fname)
{
    std::ifstream f(fname);
    if (!f) throw std::runtime_error("readGri2D: cannot open file: " + fname);

    int Nn=0, Ne=0, dim=0;
    {
        std::string line;
        std::getline(f, line);
        std::istringstream iss(line);
        iss >> Nn >> Ne >> dim;
    }
    if (dim != 2) throw std::runtime_error("readGri2D: expected dim=2");

    GriMesh2D mesh;
    mesh.V.resize(Nn, 2);
    for (int i=0; i<Nn; ++i) {
        std::string line; std::getline(f, line);
        std::istringstream iss(line);
        iss >> mesh.V(i,0) >> mesh.V(i,1);
    }

    int NB=0;
    {
        std::string line; std::getline(f, line);
        std::istringstream iss(line);
        iss >> NB;
    }

    mesh.B.resize(NB);
    for (int b=0; b<NB; ++b) {
        std::string header;
        std::getline(f, header);
        while (header.size()==0 && f) std::getline(f, header);

        std::istringstream iss(header);
        int Nb=0; int tag=0; std::string name;
        iss >> Nb >> tag >> name;

        mesh.B[b].tag = tag;
        mesh.B[b].name = name;
        mesh.B[b].edges.resize(Nb, 2);

        for (int i=0; i<Nb; ++i) {
            int a=0,c=0;
            f >> a >> c;
            mesh.B[b].edges(i,0) = a-1;
            mesh.B[b].edges(i,1) = c-1;
        }
        std::string dummy;
        std::getline(f, dummy); // consume endline
    }

    mesh.E.resize(Ne, 3);
    int Ne0 = 0;
    while (Ne0 < Ne && f) {
        std::string header;
        std::getline(f, header);
        if (!f) break;
        if (header.size()==0) continue;

        std::istringstream iss(header);
        int ne_blk=0; int tag=0; std::string etype;
        iss >> ne_blk >> tag >> etype;

        for (int i=0; i<ne_blk; ++i) {
            int a=0,b=0,c=0;
            f >> a >> b >> c;
            mesh.E(Ne0+i,0) = a-1;
            mesh.E(Ne0+i,1) = b-1;
            mesh.E(Ne0+i,2) = c-1;
        }
        std::string dummy;
        std::getline(f, dummy); // consume endline

        Ne0 += ne_blk;
    }
    if (Ne0 != Ne) throw std::runtime_error("readGri2D: element count mismatch.");

    return mesh;
}

void writeGri2D(const std::string& fname, const GriMesh2D& mesh)
{
    std::ofstream f(fname);
    if (!f) throw std::runtime_error("writeGri2D: cannot open file: " + fname);

    const int Nn = static_cast<int>(mesh.V.rows());
    const int Ne = static_cast<int>(mesh.E.rows());

    f << Nn << " " << Ne << " 2\n";
    for (int i=0; i<Nn; ++i)
        f << mesh.V(i,0) << " " << mesh.V(i,1) << "\n";

    f << mesh.B.size() << "\n";
    for (const auto& blk : mesh.B) {
        const int Nb = static_cast<int>(blk.edges.rows());
        f << Nb << " " << blk.tag << " " << blk.name << "\n";
        for (int i=0; i<Nb; ++i)
            f << (blk.edges(i,0)+1) << " " << (blk.edges(i,1)+1) << "\n";
    }

    // single triangle block
    f << Ne << " 1 TriLagrange\n";
    for (int e=0; e<Ne; ++e)
        f << (mesh.E(e,0)+1) << " " << (mesh.E(e,1)+1) << " " << (mesh.E(e,2)+1) << "\n";
}

} // namespace mesh
