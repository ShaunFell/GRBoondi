#ifndef PROCAFIELD_H_INCLUDED
#define PROCAFIELD_H_INCLUDED

/*
This class adds the simplest L2 lagrangian to the base equations of motion
*/
#include "BaseProcaField.hpp"
#include "KerrSchild.hpp"
#include "L2_simp.hpp"
#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"


class ProcaField: public BaseProcaField<KerrSchild>
{

    protected:

        template <class data_t>
        using MatterVars = typename ADMProcaVars::MatterVars<data_t>;

        template <class data_t>
        using MatterVarsD2 = typename ADMProcaVars::Diff2MatterVars<data_t>;

        template <class data_t>
        using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

        using L2_t = L2<DefaultG2>;


    public:

        struct params_t
        {
            double mass;
            double alpha2;
        };

        KerrSchild m_background;
        params_t m_params;
        L2_t m_L2;



        ProcaField(KerrSchild a_background, params_t a_params): BaseProcaField<KerrSchild>(a_background),  m_background(a_background), m_params(a_params)
        {
            //set up the L2 lagrangian

            DefaultG2::params_t G2_params{m_params.mass}; //Initialize G2 function parameters
            L2_t::params_t L2_params{m_params.alpha2}; //Initialize L2 Lagrangian parameters

            DefaultG2 m_G2(G2_params);
            this -> m_L2 = L2_t(m_G2, L2_params);
        };

        template <class data_t>
        void compute_emtensor_modification(
            emtensor_t<data_t> &base_emtensor, //pass by reference to allow modifications
            const MatterVars<data_t> &matter_vars,
            const MetricVars<data_t> &metric_vars,
            const MatterVars<Tensor<1, data_t>> &d1,
            const MatterVarsD2<Tensor<2, data_t>> &d2, //2nd derivs
            const MatterVars<data_t> &advec //value of the beta^i d_i(var) terms
        ) const override 
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
        ) const override
        {
            m_L2.matter_rhs_modification(total_rhs, vars, metric_vars, d1, d2, advec);
        };

    
};

#endif //PROCAFIELD_H_INCLUDED