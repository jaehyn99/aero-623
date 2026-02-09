#pragma once
#include "state.hpp"

class Flux {
public:
    virtual State operator()(const State& uL,
                             const State& uR,
                             double gamma) const = 0;
    virtual ~Flux() = default;
};
