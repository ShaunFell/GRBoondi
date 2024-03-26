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
#include "VarsTools.hpp"
#include "simd.hpp"

#include "DiagnosticVariables.hpp" //we need NUM_DIAGNOSTIC_VARS
#include "ADMProcaVars.hpp"

//! Does excision for fixed BG BH solutions
//! Note that it is does not using simd so one must set disable_simd()
template <class matter_t, class background_t>
class ExcisionDiagnostics
{
  protected:
      const double m_dx; //grid spacing
      const std::array<double, CH_SPACEDIM> m_center; //center of BH
        
      background_t m_background;
      double m_buffer;
      
      using Vars = typename ADMProcaVars::MatterVars<double>;

  public:

        //constructor
        ExcisionDiagnostics(background_t a_background, const double a_dx, const std::array<double, CH_SPACEDIM> a_center, double a_buffer=0.97): m_background{a_background}, m_dx{a_dx}, m_center{a_center}, m_buffer{a_buffer} {};

        void compute(const Cell<double> current_cell) const
        {
            //where are we?!
            const Coordinates<double> coords(current_cell, m_dx, m_center);

            //instance of matter variables
            Vars matter_vars;

            // check if we're in horizon region, with a buffer
            if (m_background.check_if_excised(coords, m_buffer) )
            {
                  //diagnostic variables have a corresponding enum, so we just pass the enum value to store_vars and iterate over all enums in DiagnosticVariables.hpp
                  //This allows for easy addition of new diagnostic variables
                  for ( int diag_var {0}; diag_var < NUM_DIAGNOSTIC_VARS; diag_var++)
                  {
                      current_cell.store_vars(0.0, diag_var);
                  }
            }

        }//end of method def

};

#endif /* EXCISIONDIAGNOSTICS_HPP_ */
