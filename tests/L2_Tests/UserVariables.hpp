/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "DiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    c_phi, //scalar part of proca field

    c_Avec1, //spatial part of proca field
    c_Avec2,
    c_Avec3,

    c_Evec1, //conjugate momentum of proca field
    c_Evec2,
    c_Evec3,

    c_Z, //auxiliary scalar field

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS> variable_names = {
    "phi",

    "Avec1", "Avec2", "Avec3",

    "Evec1", "Evec2", "Evec3",

    "Z"
    };

} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
