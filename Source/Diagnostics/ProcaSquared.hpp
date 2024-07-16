/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

/*
This diagnostic class calculates the square of the Proca field
*/

#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

template <class background_t> class ProcaSquared
{
  protected:
    const double m_dx;
    const std::array<double, CH_SPACEDIM> m_center;
    const background_t m_background;

    // Use the variable definition in ADMVars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    template <class data_t> using MatterVars = ADMProcaVars::MatterVars<data_t>;

  public:
    ProcaSquared(double a_dx, const std::array<double, CH_SPACEDIM> a_center,
                 const background_t a_background)
        : m_dx(a_dx), m_center{a_center}, m_background{
                                                             a_background} {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // compute background variables
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        // load variables from Chombo grid
        const auto matter_vars = current_cell.template load_vars<MatterVars>();

        // compute contravariant conformal spatial metric
        const auto gamma_UU{
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma)};

        data_t Asquared;
        Asquared = -matter_vars.phi * matter_vars.phi;
        FOR2(i, j)
        {
            Asquared +=
                gamma_UU[i][j] * matter_vars.Avec[i] * matter_vars.Avec[j];
        };

        current_cell.store_vars(Asquared, c_Asquared);
    };
};