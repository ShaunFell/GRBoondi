#ifndef PROCAFIELD_H_INCLUDED
#define PROCAFIELD_H_INCLUDED

/*
This class adds the L3 lagrangian to the base equations of motion
*/
#include "BaseProcaField.hpp"
#include "KerrSchild.hpp"
#include "L3.hpp"
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

        using L3_t = L3<DefaultG>;


    public:

        struct params_t
        {
            double mass;
            double alpha2;
            double vector_damping;
        };

        KerrSchild m_background;
        params_t m_params;
        L3_t m_L3;
        DefaultG m_G3;



        ProcaField(KerrSchild a_background, params_t a_params): BaseProcaField<KerrSchild, ProcaField>(a_background),  m_background(a_background), m_params(a_params)
        {
            //set up the L2 lagrangian

            DefaultG::params_t G_params{m_params.mass}; //Initialize G2 function parameters
            L3_t::params_t L3_params{m_params.alpha2}; //Initialize L2 Lagrangian parameters

            DefaultG a_G3(G_params);
            this -> m_L3 = L3_t(a_G3, L3_params);
            this -> m_G3 = a_G3;
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
            m_L3.compute_emtensor_modification(base_emtensor, vars, metric_vars, d1,  d2, advec);
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
            //add modifications coming from L2 lagrangian
            m_L3.matter_rhs_modification(total_rhs, matter_vars, metric_vars, d1, d2, advec);

            //add auxiliary field damping to minimize violation of gauss constraint

            Tensor<2, data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
            auto chris_phys { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };

            //compute G2 function and its derivatives
            data_t V { 0 };
            data_t dVdA { 0 };
            data_t dVddA { 0 };
            m_G3.compute_function(V, dVdA, dVddA, matter_vars, metric_vars);

            //Compute constraint algebra term
            data_t gnn{dVdA - 2.0 * dVddA * matter_vars.phi * matter_vars.phi};

            //covariatn derivative of vector field
            Tensor<2, data_t> DA;
            FOR2(i, j)
            {
                DA[i][j] = d1.Avec[j][i];
                FOR1(k) { 
                    DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; 
                }
            }

            //covariant divergence of vector field
            data_t DA_div { 0. };
            FOR2(i,j)
            {
                DA_div += gamma_UU[i][j] * DA[i][j];
            }

            //covariant divergence of Electric field
            Tensor<2,data_t> DE;
            FOR2(i, j)
            {
                DE[i][j] = d1.Evec[j][i];
                FOR1(k) { 
                    DE[i][j] += chris_phys[j][i][k] * matter_vars.Avec[k]; 
                }
            }


            //Add evolution of auxiliary scalar field
             total_rhs.Z += - metric_vars.lapse * m_params.vector_damping * matter_vars.Z + advec.Z + 2 * metric_vars.lapse *  m_params.alpha2 * dVdA * metric_vars.K * matter_vars.phi * matter_vars.phi - 2 * m_params.alpha2 * metric_vars.lapse * dVdA * matter_vars.phi * DA_div;
            FOR1(i)
            {
                total_rhs.Z +=  metric_vars.lapse * DE[i][i] + 2 * m_params.alpha2 * metric_vars.lapse * dVdA * matter_vars.Evec[i] * matter_vars.Avec[i];
                
                FOR1(j)
                {
                    total_rhs.Z += 2 * m_params.alpha2 * metric_vars.lapse * dVdA * gamma_UU[i][j] * matter_vars.Avec[i] * d1.phi[j];

            //Add damping term to electric field evolution
                    total_rhs.Evec[i] += metric_vars.lapse * gamma_UU[i][j] * d1.Z[j];

                    FOR2(k,l)
                    {
                        total_rhs.Z += - 2 * m_params.alpha2 * metric_vars.lapse * dVdA * gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.K_tensor[i][j] * matter_vars.Avec[k] * matter_vars.Avec[l];
                    }
                }


            }

            //Add damping term to scalar field evolution
            total_rhs.phi +=  - metric_vars.lapse * matter_vars.Z * m_params.mass * m_params.mass / (2 * gnn);
           
        }
    
};

#endif //PROCAFIELD_H_INCLUDED