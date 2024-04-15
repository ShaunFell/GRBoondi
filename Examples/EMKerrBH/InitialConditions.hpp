#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

/*
Example class to calculate initial data
This example sets the initial data to be a uniform magnetic field
*/

#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "Cell.hpp"
#include "KerrSchild.hpp" //background class
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

class Initial_EM_Conditions
{

  public:
    struct params_t
    {
        double init_amplitude;
    };

  protected:
    const params_t m_params;
    const KerrSchild::params_t m_Kerr_params;
    const double m_dx;

    template <class data_t>
    using MatterVars = typename ADMProcaVars::MatterVars<data_t>;

  public:
    Initial_EM_Conditions(const double a_dx, params_t a_params,
                          const KerrSchild::params_t a_Kerr_params)
        : m_params{a_params}, m_dx{a_dx}, m_Kerr_params{a_Kerr_params} {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // location of cell
        Coordinates<data_t> coords(current_cell, m_dx, m_Kerr_params.center);

        data_t radius{coords.get_radius()};
        double amp = m_params.init_amplitude;

        MatterVars<data_t> matter_vars;
        ;
        VarsTools::assign(
            matter_vars,
            0.); // assign all matter variables in this cell to zero;

        matter_vars.phi = -amp / radius; // assign scalar part to be that of
                                         // electrically charged point particle
        matter_vars.Evec[1] = amp / radius / radius;

        current_cell.store_vars(matter_vars); // push matter vars to cell
    }
};

#endif // INITIALCONDITIONS_HPP_