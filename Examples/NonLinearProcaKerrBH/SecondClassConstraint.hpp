#ifndef SECONDCLASSCONSTRAINT_HPP_INCLUDED
#define SECONDCLASSCONSTRAINT_HPP_INCLUDED

#include "ADMProcaVars.hpp"
#include "Cell.hpp"
#include "ProcaField.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

template <class background_t> class SecondClassConstraint
{
  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const std::array<double, CH_SPACEDIM> m_center;
    const background_t m_background;
    const ProcaField m_proca_field;

    // typename the matter variables
    template <class data_t>
    using MatterVars = typename ADMProcaVars::Vars<data_t>;

    // typename the metric variables
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::Vars<data_t>;

  public:
    SecondClassConstraint(double dx,
                          const std::array<double, CH_SPACEDIM> a_center,
                          const background_t a_background,
                          const ProcaField a_proca_field)
        : m_dx{dx}, m_deriv(dx), m_background(a_background), m_center(a_center),
          m_proca_field(a_proca_field) {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        MetricVars<data_t> metric_vars;

        // where are we on the grid
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // compute the metric background
        m_background.compute_metric_background(metric_vars, coords);

        // load matter variables from Chombo grid
        const MatterVars<data_t> matter_vars =
            current_cell.template load_vars<MatterVars>();

        // compute contravariant spatial metric
        const Tensor<2, data_t> gamma_UU{
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma)};

        // compute the G2 function
        data_t V{0.0};
        data_t dVdA{0.0};
        data_t dVddA{0.0};
        m_proca_field.m_G2.compute_function(V, dVdA, dVddA, matter_vars,
                                            metric_vars);

        // compute the constraint algebra term
        data_t gnn{dVdA - 2 * matter_vars.phi * matter_vars.phi * dVddA};

        // save back to chombo grid
        current_cell.store_vars(gnn, c_gnn);
    }
};

#endif // SECONDCLASSCONSTRAINT_HPP_INCLUDED