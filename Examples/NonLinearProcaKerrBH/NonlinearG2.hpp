#ifndef NONLINEARG2_HPP_INCLUDED
#define NONLINEARG2_HPP_INCLUDED

#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "TensorAlgebra.hpp"

class NonlinearG2
{
  public:
    struct params_t
    {
        double mass;
        double self_interaction;
    };

  protected:
    params_t m_params;

    template <class data_t> using MatterVars = ADMProcaVars::Vars<data_t>;

    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

  public:
    NonlinearG2(){}; // Default constructor for default initialization in
                     // matter class

    NonlinearG2(params_t a_params) : m_params{a_params} {};

    template <class data_t, template <typename> class vars_t>
    void compute_function(data_t &g_fun, data_t &g_prime, data_t &g_prime2,
                          const vars_t<data_t> &vars,
                          const MetricVars<data_t> &metric_vars) const
    {
        auto gamma_UU{TensorAlgebra::compute_inverse_sym(metric_vars.gamma)};
        data_t A_squared{-vars.phi * vars.phi};
        FOR2(i, j)
        {
            A_squared += vars.Avec[i] * vars.Avec[j] * gamma_UU[i][j];
        }

        data_t mass{m_params.mass};
        data_t lambda{m_params.self_interaction};
        data_t prefactor{0.5 * mass * mass};

        g_fun = prefactor * (A_squared + 2 * lambda * A_squared * A_squared);
        g_prime = prefactor * (1 + 4 * lambda * A_squared);
        g_prime2 = prefactor * (4 * lambda);
    }
};

#endif // NONLINEARG2_HPP_INCLUDED