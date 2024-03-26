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

    c_Z_out,

    c_Tagging_Diagnostic,

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

        "Z",

        "Tagging_Diagnostic"
        

    };
}

#endif /* DIAGNOSTICVARIABLES_HPP */