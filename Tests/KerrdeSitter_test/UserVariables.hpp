/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "EmptyDiagnosticVariables.hpp"

/// This enum gives the index of every variable stored in the grid
enum
{
    c_phi = NUM_CCZ4_VARS, // scalar part of proca field

    c_Avec1, // spatial part of proca field
    c_Avec2,
    c_Avec3,

    c_Evec1, // conjugate momentum of proca field
    c_Evec2,
    c_Evec3,

    c_Z, // auxiliary scalar field

    c_Ham,
    c_Mom1,
    c_Mom2,
    c_Mom3,

    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    user_variable_names = {"phi",

                           "Avec1", "Avec2", "Avec3",

                           "Evec1", "Evec2", "Evec3",

                           "Z",

                           "Ham",

                           "Mom1",  "Mom2",  "Mom3"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, user_variable_names);

} // namespace UserVariables

#endif /* USERVARIABLES_HPP */
