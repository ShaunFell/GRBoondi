/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONTEST_HPP_
#define EXCISIONTEST_HPP_

#include "CCZ4.hpp"
// #include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t> class ExcisionTest
{
    // Use matter_t class
    using Vars = typename matter_t::template Vars<double>;

  protected:
    const double m_dx;                              //!< The grid spacing
    const std::array<double, CH_SPACEDIM> m_center; //!< The BH center
    const background_t m_background;

  public:
    ExcisionTest(const double a_dx,
                 const std::array<double, CH_SPACEDIM> a_center,
                 background_t a_background)
        : m_dx(a_dx), m_center(a_center), m_background(a_background)
    {
    }

    void compute(Cell<double> current_cell) const
    {
        const Coordinates<double> coords(current_cell, m_dx, m_center);
        bool do_I_excise = m_background.check_if_excised(coords);

        if (do_I_excise)
        {
            // the matter rhs vars within the excision zone
            // recalculate them - for now set to zero
            Vars matter_vars;
            matter_vars.enum_mapping(
                [this, &current_cell](const int &ivar, double &var)
                { var = 0.0; });

            CCZ4::Vars<double> ccz4_vars;
            ccz4_vars.enum_mapping(
                [this, &current_cell](const int &ivar, double &var)
                { var = 0.0; });

            // assign values of rhs in output box
            // also zero the constraints
            current_cell.store_vars(matter_vars);
            current_cell.store_vars(ccz4_vars);
            current_cell.store_vars(0.0, c_Ham);
            current_cell.store_vars(0.0, c_Mom1);
            current_cell.store_vars(0.0, c_Mom2);
            current_cell.store_vars(0.0, c_Mom3);
        } // else do nothing
    }
};

#endif /* EXCISIONTEST_HPP_ */
