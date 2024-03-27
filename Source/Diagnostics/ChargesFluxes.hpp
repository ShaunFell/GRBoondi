#ifndef DIAGNOSTIC_H_INCLUDED
#define DIAGNOSTIC_H_INCLUDED

#include "CoordinateTransformations.hpp"
#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "DiagnosticVariables.hpp"

template <class matter_t, class background_t>
class ChargesFluxes
{
    protected:

        // Use the variable definition in ADMFixedBGVars
        template <class data_t> 
        using MetricVars = ADMFixedBGVars::Vars<data_t>;

        template <class data_t>
        using MatterVars = ADMProcaVars::MatterVars<data_t>;

        template <class data_t>
        using MatterDiff2Vars = ADMProcaVars::Diff2MatterVars<data_t>;

        const matter_t m_matter;
        const double m_dx;
        const std::array<double, CH_SPACEDIM> m_center;
        const FourthOrderDerivatives m_deriv;
        background_t m_background;

    public:
        ChargesFluxes(background_t a_background, double a_dx, matter_t a_matter, std::array<double, CH_SPACEDIM> a_center): m_background{a_background}, m_matter{a_matter}, m_dx{a_dx}, m_center{a_center}, m_deriv{a_dx} {};

        template <class data_t>
        void compute(Cell<data_t> current_cell) const;
};

#include "ChargesFluxes.impl.hpp"
#endif //DIAGNOSTIC_H_INCLUDED