#ifndef BASEPROCAFIELD_H_INCLUDED
#define BASEPROCAFIELD_H_INCLUDED

#include "ADMFixedBGVars.hpp"         //For metric variables
#include "ADMProcaVars.hpp"           //For matter variables
#include "DefaultBackground.hpp"      //Minkowski background as default
#include "FourthOrderDerivatives.hpp" //For calculating derivatives
#include "Tensor.hpp"                 //For performing tensorial operations
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //For user-defined variables (e.g. see EMKerrBH)
#include "simd.hpp"          //for SIMD operations

template <class background_t, class modification_t> class BaseProcaField
{
  protected:
    const background_t m_background;

    struct params_t
    {
    };

  public:
    // constructor, inputs are matter params
    BaseProcaField(background_t a_background) : m_background{a_background} {};

    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    template <class data_t> using MatterVars = ADMProcaVars::MatterVars<data_t>;

    template <class data_t>
    using Diff2MatterVars = ADMProcaVars::Diff2MatterVars<data_t>;

    /* NOTE:
    We use the method of 'Curiously Recurring Template Pattern' to allow
    arbitrary templated modifications to the theory
    https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

    This allows us to use virtual functions that can be specified by derived
    classes. e.g. see ProcaKerrBH Example
    */

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &matter_vars, // the value fo the variables
        const MetricVars<data_t> &metric_vars,
        const vars_t<Tensor<1, data_t>> &d1,       // value of 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, // 2nd derivs
        const vars_t<data_t> &advec // value of the beta^i d_i(var) terms
    ) const;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void
    matter_rhs(rhs_vars_t<data_t> &total_rhs,     // RHS terms for all vars
               const vars_t<data_t> &matter_vars, // the value fo the variables
               const MetricVars<data_t> &metric_vars,
               const vars_t<Tensor<1, data_t>> &d1,       // value of 1st derivs
               const diff2_vars_t<Tensor<2, data_t>> &d2, // 2nd derivs
               const vars_t<data_t> &advec // value of the beta^i d_i(var) terms
    ) const;
};

#include "BaseProcaField.impl.hpp"
#endif // BASEPROCAFIELD_H_INCLUDED
