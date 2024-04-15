/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

/*
enum holding diagnostic variables
*/
#ifndef DIAGNOSTICVARIABLES_HPP
#define DIAGNOSTICVARIABLES_HPP

#include <algorithm> // for std::find
#include <array>
#include <string>
#include <vector>

#include "VariableType.hpp" //GRChombo include

enum
{
    c_Asquared,

    c_rho,
    c_rhoE,
    c_rhoJ,

    c_chi,

    c_gnn,

    c_Z_out,

    c_EM_trace,
    c_EM_squared,

    c_rhoLinMom,
    c_sourceLinMom,

    // variables to be extracted on a surface
    c_Edot,
    c_Jdot,
    c_fluxLinMom,

    NUM_DIAGNOSTIC_VARS
};

namespace DiagnosticVariables
{
static const std::array<std::string, NUM_DIAGNOSTIC_VARS> variable_names = {

    "Asquared",

    "rho", "rhoE", "rhoJ",

    "chi",

    "gnn",

    "Z_out",

    "EMtrace", "EMsquared",

    "rhoLinMom", "sourceLinMom",

    // variables to be extracted on a surface
    "Edot", "Jdot", "fluxLinMom"};

template <typename T>
static const std::vector<T>
convert_pairs_to_enum(std::vector<std::pair<T, VariableType>> a_plot_pairs)
{
    std::vector<T> enum_array;
    for (int i = 0; i < a_plot_pairs.size(); i++)
    {
        enum_array.push_back(a_plot_pairs[i].first);
    }
    return enum_array;
}

static bool is_variable_to_extract(int var, std::vector<int> var_to_plot)
{
    return std::find(var_to_plot.begin(), var_to_plot.end(), var) !=
           var_to_plot.end();
}
} // namespace DiagnosticVariables

#endif /* DIAGNOSTICVARIABLES_HPP */