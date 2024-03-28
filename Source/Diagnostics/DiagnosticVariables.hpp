/*
enum holding diagnostic variables
*/
#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

enum
{
    c_Asquared,
    
    c_rho,
    c_rhoE,
    c_rhoJ,

    c_Edot,
    c_Jdot,

    c_chi,

    c_gnn,

    c_Z_out,

    c_EM_trace,
    c_EM_squared,

    c_rhoLinMom, 
    c_fluxLinMom,
    c_sourceLinMom,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
    static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {

        "Asquared",

        "rho",
        "rhoE",
        "rhoJ",

        "Edot",
        "Jdot",

        "chi",

        "gnn",

        "Z",

        "EMtrace",
        "EMsquared",

        "rhoLinMom",
        "fluxLinMom",
        "sourceLinMom"
    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */