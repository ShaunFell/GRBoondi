/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DUMBEXCISIONEVOLUTION_HPP_
#define DUMBEXCISIONEVOLUTION_HPP_

#include "ADMProcaVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"

//! Does excision for fixed BG solutions
class DumbExcisionEvolution
{
    // Use matter_t class
    using Vars = typename ADMProcaVars::MatterVars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const double m_inner_radius;
    const double m_outer_radius;
    const double m_outer_horizon;

  public:
    DumbExcisionEvolution(const double a_dx,
                          const std::array<double, CH_SPACEDIM> a_center,
                          double inner_radius, double outer_radius,
                          double a_horizon)
        : m_dx(a_dx), m_center(a_center), m_inner_radius(inner_radius),
          m_outer_radius(outer_radius), m_outer_horizon(a_horizon)
    {
    }

    void compute(const Cell<double> current_cell) const
    {
        // where are we?!
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        const double radius{coords.get_radius()};

        bool is_excised{(radius < m_inner_radius) || (radius > m_outer_radius)};
        bool is_Z_excised{radius < m_outer_horizon};

        if (is_excised)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars vars;
            VarsTools::assign(vars, 0.0);
            // assign values of rhs or vars in output box
            current_cell.store_vars(vars);
        } // else do nothing

        if (is_Z_excised)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to decay to zero
            Vars vars;
            VarsTools::assign(vars, 0.0);
            // assign values of rhs or vars in output box
            current_cell.store_vars(vars);
        }
    }
};

#endif /* DUMBEXCISIONEVOLUTION_HPP_ */
