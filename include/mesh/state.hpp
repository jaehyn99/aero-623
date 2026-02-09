#pragma once

#include <ostream>
#include <cmath>

struct State {
    double rho{0.0};
    double mom{0.0};
    double E{0.0};
};

inline State operator+(const State& a, const State& b) {
    return {
        a.rho + b.rho,
        a.mom + b.mom,
        a.E   + b.E
    };
}

inline State operator-(const State& a, const State& b) {
    return {
        a.rho - b.rho,
        a.mom - b.mom,
        a.E   - b.E
    };
}

inline State operator*(double s, const State& a) {
    return {
        s * a.rho,
        s * a.mom,
        s * a.E
    };
}

inline State operator*(const State& a, double s) {
    return s * a;
}

inline State& operator+=(State& a, const State& b) {
    a.rho += b.rho;
    a.mom += b.mom;
    a.E   += b.E;
    return a;
}

inline State& operator-=(State& a, const State& b) {
    a.rho -= b.rho;
    a.mom -= b.mom;
    a.E   -= b.E;
    return a;
}

