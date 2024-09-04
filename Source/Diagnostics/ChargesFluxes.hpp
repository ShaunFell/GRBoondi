/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/
/*
This class computes the various charges and fluxes most common to GRBoondi
simulations Current calculated quantities:
    - Conserved energy density
    - Eulerian energy density
    - Conserved angular momentum density (along z-axis)
    - Flux of conserved energy density
    - Flux of conserved angular momentum density
    - Trace of energy-momentum tensor
    - Square of energy-momentum tensor

See https://arxiv.org/pdf/2104.13420.pdf for definition of conserved quantities
*/
#ifndef DIAGNOSTIC_H_INCLUDED
#define DIAGNOSTIC_H_INCLUDED

#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "CoordinateTransformations.hpp"
#include "DiagnosticVariables.hpp"

template <class matter_t, class background_t> class ChargesFluxes
{
  protected:
    // Use the variable definition in ADMFixedBGVars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    template <class data_t> using MatterVars = ADMProcaVars::Vars<data_t>;

    template <class data_t>
    using MatterDiff2Vars = ADMProcaVars::Diff2Vars<data_t>;

    const matter_t m_matter;
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const FourthOrderDerivatives m_deriv;
    background_t m_background;

  public:
    ChargesFluxes(background_t a_background, double a_dx, matter_t a_matter,
                  std::array<double, CH_SPACEDIM> a_center)
        : m_background{a_background}, m_matter{a_matter}, m_dx{a_dx},
          m_center{a_center}, m_deriv{a_dx} {};

    template <class data_t> void compute(Cell<data_t> current_cell) const;
};

#include "ChargesFluxes.impl.hpp"
#endif // DIAGNOSTIC_H_INCLUDED