/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 *
 * Modified by GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

#ifndef EXCISIONEVOLUTION_HPP_
#define EXCISIONEVOLUTION_HPP_

#include "ADMProcaVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"

//! Does excision for fixed BG solutions
template <class matter_t, class background_t> class ExcisionEvolution
{
    // Use matter_t class
    using Vars = typename ADMProcaVars::MatterVars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const FourthOrderDerivatives m_deriv;
    const background_t m_background;
    const double m_buffer;

  public:
    ExcisionEvolution(const double a_dx,
                      const std::array<double, CH_SPACEDIM> a_center,
                      double a_buffer, background_t a_background)
        : m_dx(a_dx), m_deriv(m_dx), m_center(a_center), m_buffer(a_buffer),
          m_background(a_background)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        // where are we?!
        const Coordinates<double> coords(current_cell, m_dx, m_center);

        // check if we're inside buffered horizon region
        bool is_excised = m_background.check_if_excised(coords, m_buffer);

        std::cout << "I'm in here" << std::endl;

        if (is_excised)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars vars;
            VarsTools::assign(vars, 0.0);
            // assign values of rhs or vars in output box
            current_cell.store_vars(vars);
        } // else do nothing

        // Excise Z field right at horizon since it can derive errors
        bool is_excised_Z = m_background.check_if_excised(coords, 1.0);
        if (is_excised_Z)
        {
            current_cell.store_vars(0.0, c_Z);
        };
    }
};

#endif /* EXCISIONEVOLUTION_HPP_ */
