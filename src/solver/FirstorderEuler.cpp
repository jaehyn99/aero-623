#include "FirstorderEuler.h"
#include "mesh/TriangularMesh.h"

// Project assumption: Eigen is available, so Roe/HLLE flux headers are used directly.
#include "solver/hlle_flux.hpp"
#include "solver/roe_flux.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <iostream>
#include <limits>
#include <memory>
#include <unordered_map>
#include <string>
#include <sstream>
#include <stdexcept>

namespace {

std::string stripExtension(const std::string& path) {
    const std::size_t dot = path.find_last_of('.');
    if (dot == std::string::npos) {
        return path;
    }
    return path.substr(0, dot);
}

std::string toLower(std::string s) {
    for (char& c : s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}


EulerBoundaryConditions makeBoundaryConditions(const FirstorderEuler::SolverConfig& config) {
    EulerBoundaryConditions::Config bc;
    bc.gamma = config.gamma;
    bc.gasConstant = config.gasConstant;
    bc.alpha = config.alpha;
    bc.pout = config.pout;
    bc.rho0 = config.rho0;

    // Total-condition defaults from existing solver config.
    bc.Tt = (config.a0 * config.a0) / (config.gamma * config.gasConstant);
    bc.pt = config.rho0 * config.gasConstant * bc.Tt;
    bc.Vrot = config.a0;
    return EulerBoundaryConditions(bc);
}

double spectralRadius(const FirstorderEuler::Conserved& U,
                      const FirstorderEuler::Vec2& n,
                      double gamma);

Eigen::Vector4d toEigen(const FirstorderEuler::Conserved& u) {
    return Eigen::Vector4d(u[0], u[1], u[2], u[3]);
}

FirstorderEuler::Conserved fromEigen(const Eigen::Vector4d& u) {
    return {u(0), u(1), u(2), u(3)};
}

std::unique_ptr<Flux> makeFlux(const std::string& schemeName) {
    const std::string name = toLower(schemeName);
    if (name == "roe") {
        return std::make_unique<RoeFlux>();
    }
    if (name == "hlle") {
        return std::make_unique<HLLEFlux>();
    }
    throw std::runtime_error("Unknown flux scheme: " + schemeName + " (expected 'roe' or 'hlle').");
}

double spectralRadius(const FirstorderEuler::Conserved& U,
                      const FirstorderEuler::Vec2& n,
                      double gamma) {
    const double rho = std::max(1e-14, U[0]);
    const double u = U[1] / rho;
    const double v = U[2] / rho;
    const double p = std::max(1e-14, (gamma - 1.0) * (U[3] - 0.5 * rho * (u * u + v * v)));
    const double c = std::sqrt(gamma * p / rho);
    const double un = u * n[0] + v * n[1];
    return std::abs(un) + c;
}

} // namespace

FirstorderEuler::MeshInputs FirstorderEuler::MeshInputs::fromMeshFile(const std::string& meshFile) {
    return fromPrefix(stripExtension(meshFile));
}

FirstorderEuler::MeshInputs FirstorderEuler::MeshInputs::fromPrefix(const std::string& prefix) {
    MeshInputs inputs;
    inputs.meshFile = prefix + ".gri";
    inputs.b2eFile = prefix + "B2E.txt";
    inputs.bnFile = prefix + "Bn.txt";
    inputs.i2eFile = prefix + "I2E.txt";
    inputs.inFile = prefix + "In.txt";
    inputs.periodicEdgesFile = prefix + "periodicEdges.txt";
    return inputs;
}

FirstorderEuler::FirstorderEuler(MeshInputs inputs, SolverConfig config)
    : inputs_(std::move(inputs)), config_(std::move(config)) {}

void FirstorderEuler::loadInputs() {
    readMeshAndConnectivity();
    if (config_.validateMeshOnLoad) {
        validateLoadedArrays();
    }

    U_.assign(elements_.size(), Conserved{0.0, 0.0, 0.0, 0.0});
    residual_.assign(elements_.size(), Conserved{0.0, 0.0, 0.0, 0.0});
}

void FirstorderEuler::readMeshAndConnectivity() {
    // Simpler flow: parse .gri once with TriangularMesh and derive element/face arrays from it.
    triMesh_ = std::make_unique<TriangularMesh>(inputs_.meshFile);

    nodes_.clear();
    nodes_.reserve(static_cast<std::size_t>(triMesh_->numNodes()));
    for (std::size_t i = 0; i < static_cast<std::size_t>(triMesh_->numNodes()); ++i) {
        const auto& p = triMesh_->node(i);
        nodes_.push_back({p[0], p[1]});
    }

    elements_.clear();
    elements_.reserve(static_cast<std::size_t>(triMesh_->numElems()));
    for (std::size_t i = 0; i < static_cast<std::size_t>(triMesh_->numElems()); ++i) {
        const auto& e = triMesh_->elem(i);
        elements_.push_back({static_cast<std::size_t>(e._pointID[0]),
                             static_cast<std::size_t>(e._pointID[1]),
                             static_cast<std::size_t>(e._pointID[2])});
    }

    area_.clear();
    area_.reserve(static_cast<std::size_t>(triMesh_->numElems()));
    for (std::size_t i = 0; i < static_cast<std::size_t>(triMesh_->numElems()); ++i) {
        const double ai = triMesh_->area(i);
        if (ai <= 0.0) {
            throw std::runtime_error("TriangularMesh area is non-positive.");
        }
        area_.push_back(ai);
    }

    buildFacesFromTriangularMesh();
    computePerimeterFromMesh();
}

void FirstorderEuler::buildFacesFromTriangularMesh() {
    if (!triMesh_) {
        throw std::runtime_error("TriangularMesh must be initialized before building face data.");
    }

    interiorFaces_.clear();
    boundaryFaces_.clear();
    periodicEdges_.clear();

    std::unordered_map<std::string, std::size_t> groupFromTitle;

    for (std::size_t faceID = 0; faceID < static_cast<std::size_t>(triMesh_->numFaces()); ++faceID) {
        const auto& face = triMesh_->face(faceID);
        const bool periodic = face._periodicFaceID != -1;
        const bool boundary = face.isBoundaryFace() && !periodic;

        const std::size_t elemL = static_cast<std::size_t>(face._elemID[0]);
        const std::size_t localFaceL = static_cast<std::size_t>(std::find(
            triMesh_->elem(elemL)._faceID.begin(),
            triMesh_->elem(elemL)._faceID.end(),
            static_cast<int>(faceID)) - triMesh_->elem(elemL)._faceID.begin());

        if (boundary) {
            const auto [it, inserted] = groupFromTitle.try_emplace(face._title, groupFromTitle.size());
            (void)inserted;

            BoundaryFace bf;
            bf.elem = elemL;
            bf.localFace = localFaceL;
            bf.boundaryGroup = it->second;
            bf.boundaryTitle = face._title;
            bf.normal = {triMesh_->normal(elemL, localFaceL)[0], triMesh_->normal(elemL, localFaceL)[1]};
            const auto& p0 = triMesh_->node(static_cast<std::size_t>(face._pointID[0]));
            const auto& p1 = triMesh_->node(static_cast<std::size_t>(face._pointID[1]));
            bf.center = {(p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5};
            bf.length = face._length;
            boundaryFaces_.push_back(bf);
            continue;
        }

        std::size_t elemR = 0;
        std::size_t localFaceR = 0;
        if (periodic) {
            elemR = static_cast<std::size_t>(face._periodicElemID);
            localFaceR = static_cast<std::size_t>(std::find(
                triMesh_->elem(elemR)._faceID.begin(),
                triMesh_->elem(elemR)._faceID.end(),
                face._periodicFaceID) - triMesh_->elem(elemR)._faceID.begin());
        } else {
            elemR = static_cast<std::size_t>(face._elemID[1]);
            localFaceR = static_cast<std::size_t>(std::find(
                triMesh_->elem(elemR)._faceID.begin(),
                triMesh_->elem(elemR)._faceID.end(),
                static_cast<int>(faceID)) - triMesh_->elem(elemR)._faceID.begin());
        }

        if (elemR < elemL) {
            continue;
        }

        InteriorFace inf;
        inf.elemL = elemL;
        inf.faceL = localFaceL;
        inf.elemR = elemR;
        inf.faceR = localFaceR;
        inf.normal = {triMesh_->normal(elemL, localFaceL)[0], triMesh_->normal(elemL, localFaceL)[1]};
        const auto& p0 = triMesh_->node(static_cast<std::size_t>(face._pointID[0]));
        const auto& p1 = triMesh_->node(static_cast<std::size_t>(face._pointID[1]));
        inf.center = {(p0[0] + p1[0]) * 0.5, (p0[1] + p1[1]) * 0.5};
        inf.length = face._length;
        interiorFaces_.push_back(inf);

        if (periodic) {
            periodicEdges_.push_back({elemL, elemR});
        }
    }
}

void FirstorderEuler::computePerimeterFromMesh() {
    if (!triMesh_) {
        throw std::runtime_error("TriangularMesh must be initialized before computing perimeter.");
    }

    perimeter_.assign(elements_.size(), 0.0);
    for (std::size_t elemID = 0; elemID < static_cast<std::size_t>(triMesh_->numElems()); ++elemID) {
        const auto& elem = triMesh_->elem(elemID);
        perimeter_[elemID] = triMesh_->length(static_cast<std::size_t>(elem._faceID[0]))
            + triMesh_->length(static_cast<std::size_t>(elem._faceID[1]))
            + triMesh_->length(static_cast<std::size_t>(elem._faceID[2]));
    }
}

void FirstorderEuler::validateLoadedArrays() const {
    if (nodes_.empty() || elements_.empty()) {
        throw std::runtime_error("Mesh must contain nodes and elements.");
    }
    if (area_.size() != elements_.size()) {
        throw std::runtime_error("TriangularMesh area count must equal number of elements.");
    }

    for (const auto& f : interiorFaces_) {
        if (f.elemL >= elements_.size() || f.elemR >= elements_.size()) {
            throw std::runtime_error("I2E element index out of range.");
        }
        if (f.faceL > 2 || f.faceR > 2) {
            throw std::runtime_error("I2E local face index must be in {1,2,3} (0-based {0,1,2}).");
        }
    }

    for (const auto& f : boundaryFaces_) {
        if (f.elem >= elements_.size()) {
            throw std::runtime_error("B2E element index out of range.");
        }
        if (f.localFace > 2) {
            throw std::runtime_error("B2E local face index must be in {1,2,3} (0-based {0,1,2}).");
        }
    }
}

void FirstorderEuler::initUniformState() {
    // Physical meaning:
    // Start from inlet stagnation conditions (rho0, a0) plus a guessed Mach number.
    // Convert them to one thermodynamically consistent static state (rho, u, v, p),
    // then initialize every cell with the same conservative vector.
    // This gives a stable and simple first iterate before edge flux updates begin.
    const double M = config_.initialMach;
    const double g = config_.gamma;
    const double R = config_.gasConstant;
    const double p0 = config_.rho0 * config_.a0 * config_.a0 / g;

    const double T0 = (config_.a0 * config_.a0) / (g * R);
    const double T = T0 / (1.0 + 0.5 * (g - 1.0) * M * M);
    const double p = p0 * std::pow(T / T0, g / (g - 1.0));

    const double u = M * config_.a0 * std::cos(config_.alpha);
    const double v = M * config_.a0 * std::sin(config_.alpha);

    const double rho = p / (R * T);
    const double rhoE = p / (g - 1.0) + 0.5 * rho * (u * u + v * v);

    for (auto& Ui : U_) {
        Ui[0] = rho;
        Ui[1] = rho * u;
        Ui[2] = rho * v;
        Ui[3] = rhoE;
    }
}

void FirstorderEuler::advanceToConvergedOrFinalTime() {
    advance(true);
}

void FirstorderEuler::runSteadyGlobal() {
    config_.localTimeStepping = false;
    resetMarchState();
    advance(false);
}

void FirstorderEuler::runSteadyLocal() {
    config_.localTimeStepping = true;
    resetMarchState();
    advance(false);
}

void FirstorderEuler::runUnsteadyGlobal() {
    config_.localTimeStepping = false;
    resetMarchState();
    advance(true);
}

void FirstorderEuler::resetMarchState() {
    time_ = 0.0;
    iteration_ = 0;
}

void FirstorderEuler::advance(bool stopByTime) {
    std::cout << "    it           t          dt        ||R||2\n";

    while (iteration_ < config_.maxIterations && (!stopByTime || time_ < config_.finalTime)) {
        auto edges = computeEdgeFluxesAndWaveSpeeds();
        assembleResidualFromEdgeFluxes(edges);

        double dtUsed = 0.0;
        if (config_.localTimeStepping) {
            auto dtLocal = computeLocalDtFromWaveSpeeds(edges);
            dtUsed = dtLocal.empty() ? 0.0 : *std::min_element(dtLocal.begin(), dtLocal.end());
            updateStateLocalDt(dtLocal);
        } else {
            dtUsed = computeGlobalDtFromWaveSpeeds(edges);
            updateStateGlobalDt(dtUsed);
        }

        const double normR = l2Norm(residual_);
        if (iteration_ % 10 == 0 || (stopByTime && time_ >= config_.finalTime)) {
            std::cout << iteration_ << " " << time_ << " " << dtUsed << " " << normR << "\n";
        }

        ++iteration_;
        time_ += dtUsed;

        if (normR < config_.residualTolerance) {
            break;
        }
    }
}

std::vector<FirstorderEuler::EdgeFluxContribution> FirstorderEuler::computeEdgeFluxesAndWaveSpeeds() const {
    std::vector<EdgeFluxContribution> edges;
    edges.reserve(interiorFaces_.size() + boundaryFaces_.size());

    const auto flux = makeFlux(config_.fluxScheme);
    const EulerBoundaryConditions bcModel = makeBoundaryConditions(config_);

    for (const auto& f : interiorFaces_) {
        const Conserved& UL = U_.at(f.elemL);
        const Conserved& UR = U_.at(f.elemR);

        const Eigen::Vector2d n(f.normal[0], f.normal[1]);
        const Conserved F = fromEigen((*flux)(toEigen(UL), toEigen(UR), config_.gamma, n));

        EdgeFluxContribution edge;
        edge.ownerElem = f.elemL;
        edge.neighborElem = f.elemR;
        edge.normal = f.normal;
        edge.edgeLength = f.length;
        edge.flux = F;
        edge.spectralRadius = std::max(spectralRadius(UL, f.normal, config_.gamma),
                                       spectralRadius(UR, f.normal, config_.gamma));
        edges.push_back(edge);
    }

    // Boundary condition imposition happens here via ghost-state construction.
    for (const auto& f : boundaryFaces_) {
        const Conserved& UL = U_.at(f.elem);
        EulerBoundaryConditions::Type kind = EulerBoundaryConditions::typeFromCurveTitle(f.boundaryTitle);

        // For unsteady inflow support, switch inflow types by runtime mode convention.
        if (kind == EulerBoundaryConditions::Type::InflowSteady
            && !config_.localTimeStepping
            && config_.finalTime < 1e11) {
            kind = EulerBoundaryConditions::Type::InflowUnsteady;
        }

        EulerBoundaryConditions::Context ctx;
        ctx.time = time_;
        ctx.faceCenter = f.center;
        const Conserved UR = bcModel.boundaryState(kind, UL, f.normal, ctx);

        const Eigen::Vector2d n(f.normal[0], f.normal[1]);
        const Conserved F = fromEigen((*flux)(toEigen(UL), toEigen(UR), config_.gamma, n));

        EdgeFluxContribution edge;
        edge.ownerElem = f.elem;
        edge.neighborElem = f.elem;
        edge.normal = f.normal;
        edge.edgeLength = f.length;
        edge.flux = F;
        edge.spectralRadius = std::max(spectralRadius(UL, f.normal, config_.gamma),
                                       spectralRadius(UR, f.normal, config_.gamma));
        edges.push_back(edge);
    }

    return edges;
}

void FirstorderEuler::assembleResidualFromEdgeFluxes(const std::vector<EdgeFluxContribution>& edges) {
    for (auto& R : residual_) {
        for (double& x : R) {
            x = 0.0;
        }
    }

    for (const auto& edge : edges) {
        auto& RL = residual_.at(edge.ownerElem);
        for (std::size_t k = 0; k < RL.size(); ++k) {
            RL[k] += edge.edgeLength * edge.flux[k];
        }

        if (edge.neighborElem != edge.ownerElem) {
            auto& RR = residual_.at(edge.neighborElem);
            for (std::size_t k = 0; k < RR.size(); ++k) {
                RR[k] -= edge.edgeLength * edge.flux[k];
            }
        }
    }
}

double FirstorderEuler::computeGlobalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const {
    if (U_.empty()) {
        return 0.0;
    }

    if (edges.empty()) {
        const double minArea = *std::min_element(area_.begin(), area_.end());
        return std::max(1e-12, config_.cfl * minArea);
    }

    double sumS = 0.0;
    for (const auto& edge : edges) {
        sumS += std::abs(edge.spectralRadius) * edge.edgeLength;
    }

    const double avgS = sumS / static_cast<double>(edges.size());
    const double minDiameter = [&]() {
        double d = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < area_.size(); ++i) {
            d = std::min(d, 2.0 * area_[i] / std::max(1e-12, perimeter_[i]));
        }
        return d;
    }();

    return config_.cfl * minDiameter / std::max(1e-12, avgS);
}

std::vector<double> FirstorderEuler::computeLocalDtFromWaveSpeeds(const std::vector<EdgeFluxContribution>& edges) const {
    std::vector<double> dt(area_.size(), 1e-3);
    if (edges.empty()) {
        return dt;
    }

    std::vector<double> edgeSpeedWeighted(area_.size(), 0.0);
    for (const auto& edge : edges) {
        edgeSpeedWeighted[edge.ownerElem] += std::abs(edge.spectralRadius) * edge.edgeLength;
        if (edge.neighborElem != edge.ownerElem) {
            edgeSpeedWeighted[edge.neighborElem] += std::abs(edge.spectralRadius) * edge.edgeLength;
        }
    }

    for (std::size_t i = 0; i < area_.size(); ++i) {
        dt[i] = 2.0 * area_[i] * config_.cfl / std::max(1e-12, edgeSpeedWeighted[i]);
    }
    return dt;
}

void FirstorderEuler::updateStateGlobalDt(double dt) {
    for (std::size_t i = 0; i < U_.size(); ++i) {
        const double scale = dt / area_[i];
        for (std::size_t k = 0; k < U_[i].size(); ++k) {
            U_[i][k] -= scale * residual_[i][k];
        }
    }
}

void FirstorderEuler::updateStateLocalDt(const std::vector<double>& dtLocal) {
    for (std::size_t i = 0; i < U_.size(); ++i) {
        const double scale = dtLocal[i] / area_[i];
        for (std::size_t k = 0; k < U_[i].size(); ++k) {
            U_[i][k] -= scale * residual_[i][k];
        }
    }
}

double FirstorderEuler::l2Norm(const std::vector<Conserved>& values) {
    double acc = 0.0;
    for (const auto& vi : values) {
        for (double x : vi) {
            acc += x * x;
        }
    }
    return std::sqrt(acc);
}
