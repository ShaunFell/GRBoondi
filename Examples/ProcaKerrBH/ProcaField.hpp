#ifndef PROCAFIELD_H_INCLUDED
#define PROCAFIELD_H_INCLUDED

/*
This class adds the simplest L2 lagrangian to the base equations of motion
*/
#include "BaseProcaField.hpp"
#include "KerrSchild.hpp"
#include "L2_simp.hpp"
#include "DefaultG.hpp"
#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"

// Note: base class BaseProcaField uses CRTP, so pass ProcaField itself as template argument
class ProcaField: public BaseProcaField<KerrSchild, ProcaField>
{

    protected:

        template <class data_t>
        using MatterVars = typename ADMProcaVars::MatterVars<data_t>;

        template <class data_t>
        using MatterVarsD2 = typename ADMProcaVars::Diff2MatterVars<data_t>;

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
        L2_t m_L2;
        DefaultG m_G2;



        ProcaField(KerrSchild a_background, params_t a_params): BaseProcaField<KerrSchild, ProcaField>(a_background),  m_background(a_background), m_params(a_params)
        {
            //set up the L2 lagrangian

            DefaultG::params_t G2_params{m_params.mass}; //Initialize G2 function parameters
            L2_t::params_t L2_params{m_params.alpha2}; //Initialize L2 Lagrangian parameters

            DefaultG a_G2(G2_params);
            this -> m_L2 = L2_t(a_G2, L2_params);
            this -> m_G2 = a_G2;
        };

        template <class data_t>
        void compute_emtensor_modification(
            emtensor_t<data_t> &base_emtensor, //pass by reference to allow modifications
            const MatterVars<data_t> &matter_vars,
            const MetricVars<data_t> &metric_vars,
            const MatterVars<Tensor<1, data_t>> &d1,
            const MatterVarsD2<Tensor<2, data_t>> &d2, //2nd derivs
            const MatterVars<data_t> &advec //value of the beta^i d_i(var) terms
        ) const 
        {
            m_L2.compute_emtensor_modification(base_emtensor, matter_vars, metric_vars, d1,  d2, advec);
        }; 

        template <class data_t,  template <typename> class rhs_vars_t>
        void matter_rhs_modification(
            rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
            const MatterVars<data_t> &vars, //the value fo the variables
            const MetricVars<data_t> &metric_vars,
            const MatterVars<Tensor<1, data_t>> &d1, //value of 1st derivs
            const MatterVarsD2<Tensor<2, data_t>> &d2, //2nd derivs
            const MatterVars<data_t> &advec //value of the beta^i d_i(var) terms
        ) const
        {
            //add modification due to L2 lagrangian
            m_L2.matter_rhs_modification(total_rhs, vars, metric_vars, d1, d2, advec);

            //add auxiliary field damping to minimize violation of gauss constraint

            Tensor<2, data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
            auto chris_phys { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };

            data_t g_func { 0 };
            data_t g_func_prime { 0 };
            data_t g_func_prime2 { 0 };
            m_G2.compute_function(g_func, g_func_prime, g_func_prime2, vars, metric_vars, d1, d2);

            data_t gnn {  m_params.alpha2 * (g_func_prime - 2 * vars.phi * vars.phi * g_func_prime2) };

            //Modify scalar part
            total_rhs.phi -= metric_vars.lapse * vars.Z * m_params.mass * m_params.mass / (2 * gnn);

            //modify electric field part
            FOR2(i,j)
            {
                total_rhs.Evec[i] += metric_vars.lapse * gamma_UU[i][j] * d1.Z[j];
            }

            //add evolution for auxiliary Z field
            total_rhs.Z = 2 * m_params.alpha2 * metric_vars.lapse * g_func_prime * vars.phi - m_params.vector_damping * metric_vars.lapse * vars.Z + advec.Z;

            FOR1(i)
            {
                total_rhs.Z += metric_vars.lapse * d1.Evec[i][i];
                FOR1(j)
                {
                    total_rhs.Z += metric_vars.lapse * chris_phys[i][i][j] * vars.Evec[j];
                }
            }
        };

    
};

#endif //PROCAFIELD_H_INCLUDED