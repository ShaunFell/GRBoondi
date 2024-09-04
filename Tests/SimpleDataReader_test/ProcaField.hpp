#ifndef PROCAFIELD_H_INCLUDED
#define PROCAFIELD_H_INCLUDED

/*
This class adds the simplest L2 lagrangian to the base equations of motion
*/
#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "DefaultG.hpp"
#include "KerrSchild.hpp"
#include "L2_simp.hpp"

// Note: base class BaseProcaField uses CRTP, so pass ProcaField itself as
// template argument
class ProcaField : public BaseProcaField<KerrSchild, ProcaField>
{

  protected:
    template <class data_t>
    using Vars = typename ADMProcaVars::Vars<data_t>;

    template <class data_t>
    using VarsD2 = typename ADMProcaVars::Diff2Vars<data_t>;

    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

    using L2_t = L2<DefaultG>;

  public:
    struct params_t
    {
        double mass;
        double alpha2;
        double vector_damping;
    };

    KerrSchild m_background;
    params_t m_params;

    ProcaField(KerrSchild a_background, params_t a_params)
        : BaseProcaField<KerrSchild, ProcaField>(a_background),
          m_background(a_background), m_params(a_params){};

    template <class data_t>
    void compute_emtensor_modification(
        emtensor_t<data_t>
            &base_emtensor, // pass by reference to allow modifications
        const Vars<data_t> &vars, const MetricVars<data_t> &metric_vars,
        const Vars<Tensor<1, data_t>> &d1,
        const VarsD2<Tensor<2, data_t>> &d2, // 2nd derivs
        const Vars<data_t> &advec // value of the beta^i d_i(var) terms
    ) const {};

    template <class data_t, template <typename> class rhs_vars_t>
    void matter_rhs_modification(
        rhs_vars_t<data_t> &total_rhs,         // RHS terms for all vars
        const Vars<data_t> &matter_vars, // the value fo the variables
        const MetricVars<data_t> &metric_vars,
        const Vars<Tensor<1, data_t>> &d1,   // value of 1st derivs
        const VarsD2<Tensor<2, data_t>> &d2, // 2nd derivs
        const Vars<data_t> &advec // value of the beta^i d_i(var) terms
    ) const {};
};

#endif // PROCAFIELD_H_INCLUDED