#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace mconst{
    // Mathematical constants
    constexpr double pi = 3.14159265358979323846;
}

namespace pconst{
    // Physical constants
    constexpr double c = 299792458; // Speed of light in vacuum [m/s]
    constexpr double e = 1.602176634e-19; // Elementary charge [C]
    constexpr double h = 6.62607015e-34; // Planck constant [J s]
    constexpr double k_B = 1.380649e-23; // Boltzmann constant [J/kg]
    constexpr double m_e = 9.10938371393e-31; // Electron mass [kg]
    constexpr double m_n = 1.67492750056e-27; // Neutron mass [kg]
    constexpr double R = 8.31446261815324; // Universal gas constant [J/mol K]
}

#endif