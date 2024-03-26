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
            const MatterVars<data_t> &vars,
            const MetricVars<data_t> &metric_vars,
            const MatterVars<Tensor<1, data_t>> &d1,
            const MatterVarsD2<Tensor<2, data_t>> &d2, //2nd derivs
            const MatterVars<data_t> &advec //value of the beta^i d_i(var) terms
        ) const 
        {
            m_L2.compute_emtensor_modification(base_emtensor, vars, metric_vars, d1,  d2, advec);
        }; 

        template <class data_t,  template <typename> class rhs_vars_t>
        void matter_rhs_modification(
            rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
            const MatterVars<data_t> &matter_vars, //the value fo the variables
            const MetricVars<data_t> &metric_vars,
            const MatterVars<Tensor<1, data_t>> &d1, //value of 1st derivs
            const MatterVarsD2<Tensor<2, data_t>> &d2, //2nd derivs
            const MatterVars<data_t> &advec //value of the beta^i d_i(var) terms
        ) const
        {
            /*
            //add modification due to L2 lagrangian
            m_L2.matter_rhs_modification(total_rhs, matter_vars, metric_vars, d1, d2, advec);

            //add auxiliary field damping to minimize violation of gauss constraint

            Tensor<2, data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
            auto chris_phys { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };

            data_t g_func { 0 };
            data_t g_func_prime { 0 };
            data_t g_func_prime2 { 0 };
            m_G2.compute_function(g_func, g_func_prime, g_func_prime2, matter_vars, metric_vars, d1, d2);

            data_t gnn {  m_params.alpha2 * (g_func_prime - 2 * matter_vars.phi * matter_vars.phi * g_func_prime2) };
            data_t mass{m_params.mass};


            // covariant derivative of spatial part of Proca field
            Tensor<2, data_t> DA;
            FOR2(i, j)
            {
                DA[i][j] = d1.Avec[j][i];
                FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
            }

           FOR1(i)
            {
                total_rhs.Avec[i] =
                    -metric_vars.lapse * d1.phi[i] - matter_vars.phi * metric_vars.d1_lapse[i] + advec.Avec[i];

                FOR1(j)
                {
                    total_rhs.Avec[i] += -metric_vars.lapse * metric_vars.gamma[i][j] * matter_vars.Evec[j] +
                                        matter_vars.Avec[j] * metric_vars.d1_shift[j][i];
                };
            };

            //modify electric field part
            FOR1(i)
            {
                total_rhs.Evec[i] = metric_vars.lapse * metric_vars.K * matter_vars.Evec[i] + advec.Evec[i];

                FOR1(j)
                {
                    total_rhs.Evec[i] += metric_vars.lapse * gamma_UU[i][j] * d1.Z[j] + m_params.alpha2 * ( 2 * metric_vars.lapse * g_func_prime * gamma_UU[i][j] * matter_vars.Avec[j]  ) - matter_vars.Evec[j] * metric_vars.d1_shift[i][j];
                }

                FOR3(j, k, l)
                {
                    total_rhs.Evec[i] += gamma_UU[j][k] * gamma_UU[i][l] * (
                                            metric_vars.d1_lapse[j] * ( d1.Avec[k][l] - d1.Avec[l][k] ) +
                                            metric_vars.lapse * ( d2.Avec[k][l][j] - d2.Avec[l][k][j] )
                                        );

                    FOR1(m)
                    {
                        total_rhs.Evec[i] +=
                            -metric_vars.lapse * gamma_UU[j][k] * gamma_UU[i][l] *
                            (chris_phys[m][j][l] * (d1.Avec[k][m] - d1.Avec[m][k]) +
                            chris_phys[m][j][k] * (d1.Avec[m][l] - d1.Avec[l][m]));
                    };
                };
            }; 

            //add evolution for auxiliary Z field
            total_rhs.Z += 2 * metric_vars.lapse * g_func_prime * matter_vars.phi - m_params.vector_damping * metric_vars.lapse * matter_vars.Z + advec.Z;
            FOR1(i)
            {
                total_rhs.Z += metric_vars.lapse * d1.Evec[i][i];
                FOR1(j)
                {
                    total_rhs.Z += metric_vars.lapse * chris_phys[i][i][j] * matter_vars.Evec[j];
                }
            } 
        
            total_rhs.phi = - metric_vars.lapse * matter_vars.Z * m_params.mass * m_params.mass / (2 * gnn) +
                                m_params.alpha2 * metric_vars.lapse * g_func_prime * matter_vars.phi * metric_vars.K / (gnn) + advec.phi;
            FOR1(i)
            {
                total_rhs.phi += 2 * metric_vars.lapse * m_params.alpha2 * g_func_prime2 * matter_vars.phi * matter_vars.Avec[i] *
                                matter_vars.Evec[i] / gnn;

                FOR1(j)
                {
                    total_rhs.phi +=
                        gamma_UU[i][j] * (-metric_vars.lapse * m_params.alpha2 * g_func_prime / gnn * DA[i][j] -
                                        matter_vars.Avec[i] * metric_vars.d1_lapse[j] +
                                        2 * metric_vars.lapse * m_params.alpha2 * g_func_prime2 / gnn * 2 * matter_vars.phi *
                                            matter_vars.Avec[i] * d1.phi[j]);

                    FOR2(k, l)
                    {
                        total_rhs.phi -=
                            gamma_UU[i][k] * gamma_UU[j][l] * m_params.alpha2 * 
                            (2 * metric_vars.lapse * g_func_prime2 / gnn * matter_vars.phi * matter_vars.Avec[i] *
                                matter_vars.Avec[j] * metric_vars.K_tensor[k][l] +
                            2 * metric_vars.lapse * g_func_prime2 / gnn * matter_vars.Avec[i] *
                                matter_vars.Avec[j] * DA[k][l]);
                    }
                }
            };

            */

        }
    
};

#endif //PROCAFIELD_H_INCLUDED