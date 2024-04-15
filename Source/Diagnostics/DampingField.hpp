/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

/*
  This class computes the value of the auxiliary Z field and assigns it to the
  diagnostic variable
*/

#ifndef DAMPINGFIELD_HPP_INCLUDED
#define DAMPINGFIELD_HPP_INCLUDED

#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

class DampingFieldDiagnostic
{
  protected:
    template <class data_t>
    using MatterVars =
        typename ADMProcaVars::MatterVars<data_t>; // type alias the matter
                                                   // variables

  public:
    DampingFieldDiagnostic(){}; // explicit default constructor

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // Load the auxiliary Z field and save it to diagnostic variable
        // c_Zvec_out. The auxiliary Z field is equal to minus the Gauss
        // constraint. so instead of calculating the Gauss constraint directly,
        // we just have to load the Z field

        MatterVars<data_t> matter_vars{current_cell.template load_vars<
            MatterVars>()}; // load the matter variables from the Chombo grid
        current_cell.store_vars(matter_vars.Z, c_Z_out);
    };
};

#endif /* DAMPINGFIELD_HPP_INCLUDED */