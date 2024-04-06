#ifndef DEFAULT_G_HPP
#define DEFAULT_G_HPP

#include "ADMProcaVars.hpp"
#include "ADMFixedBGVars.hpp"
#include "TensorAlgebra.hpp"

class DefaultG
{
    public:
        struct params_t
        {
            double mass;
        };

    
    protected:
        params_t m_params;

        template <class data_t>
        using MatterVars = ADMProcaVars::MatterVars<data_t>;

        template <class data_t>
        using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;
    
    public:
        DefaultG(){}; //Default constructor for default initialization in matter class

        DefaultG(params_t a_params): m_params{a_params} {};

        template <class data_t, template <typename> class vars_t>
        void compute_function(data_t &g_fun, data_t &g_prime, data_t &g_prime2, const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars) const
        {
            auto gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
            data_t A_squared { - vars.phi * vars.phi };
            FOR2(i,j)
            {
                A_squared += vars.Avec[i] * vars.Avec[j] * gamma_UU[i][j];
            }

            g_fun = 0.5 * m_params.mass * m_params.mass * A_squared;
            g_prime = 0.5 * m_params.mass * m_params.mass; 
            g_prime2 = 0.;
        }

};

#endif //DEFAULT_G_HPP