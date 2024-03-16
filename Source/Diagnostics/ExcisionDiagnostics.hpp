/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EXCISIONDIAGNOSTICS_HPP_
#define EXCISIONDIAGNOSTICS_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

#include "DiagnosticVariables.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t>
class ExcisionDiagnostics
{
  protected:
        const double m_dx; //grid spacing
        const std::array<double, CH_SPACEDIM> m_center; //center of BH
        
        background_t m_background;

      using Vars = typename matter_t::template Vars<double>;

  public:

        //constructor
        ExcisionDiagnostics(background_t a_background, const double a_dx, const std::array<double, CH_SPACEDIM> a_center): m_background{a_background}, m_dx{a_dx}, m_center{a_center}
        {
        };

        void compute(const Cell<double> current_cell) const
        {
            const Coordinates<double> coords(current_cell, m_dx, m_center);
            Vars matter_vars;

            if (m_background.check_if_excised(coords) )
            {
                  VarsTools::assign(matter_vars, 0.0);

                  //assign values of variables to cell
                  current_cell.store_vars(matter_vars);
            }

        }//end of method def

};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */
