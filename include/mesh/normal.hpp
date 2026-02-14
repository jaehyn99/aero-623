#pragma once
#include "state.hpp"

class Normal {
public:
    virtual State operator()(const State& uL,
                             const State& uR,
                             double gamma, const Normal& n) const = 0;
    virtual ~Flux() = default;
};
