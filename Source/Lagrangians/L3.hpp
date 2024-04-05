#ifndef L3_SIMP_H_INCLUDED
#define L3_SIMP_H_INCLUDED


/*
This file calculates to the modifications to the equations of motion and energy momentum tensor coming from the third generalized Proca lagrangian

          -1      2          
L   =  ── ⋅ F  + α  ⋅ L 
 gp     4               3    3

                            ^^ We calculate this term


IMPORTANT: This code is still untested. Use at your own risk!
                       
*/



#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "TensorAlgebra.hpp"
#include "DefaultG.hpp"

/*
This class calculates to the modifications to the equations of motion and energy momentum tensor coming from the second generalized Proca lagrangian, 
which is just an arbitrary function.


NOTE: This assumes the theory is L = -1/4 F^2 + alpha3 * L_3 with no other modifications due to higher lagrangians
*/



template <class G3 = DefaultG>
class L3
{
    public:
        struct params_t
        {
            double alpha3;
        };

    protected:

        G3 m_g3_function;
        params_t m_params;

        template <class data_t>
        using MatterVars = ADMProcaVars::MatterVars<data_t>;

        template <class data_t>
        using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

    public:

        L3() {}; //Default constructor for default initialization in matter class 

        L3(G3 a_G3_function, params_t a_params) : m_g3_function(a_G3_function), m_params{a_params} {};

    
        template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
        void compute_emtensor_modification(
            emtensor_t<data_t> &base_emtensor,
            const vars_t<data_t> &matter_vars,
            const MetricVars<data_t> &metric_vars,
            const vars_t<Tensor<1,data_t>> &d1,
            const diff2_vars_t<Tensor<2,data_t>> &d2,
            const vars_t<data_t> &advec
        ) const;


        template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
        void matter_rhs_modification(
            rhs_vars_t<data_t> &total_rhs,
            const vars_t<data_t> &vars,
            const MetricVars<data_t> &metric_vars,
            const vars_t<Tensor<1,data_t>> &d1,
            const diff2_vars_t<Tensor<2,data_t>> &d2,
            const vars_t<data_t> &advec
        ) const ;

        template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
        void compute_phi_dot(
            data_t &phi_dot,
            const vars_t<data_t> &vars,
            const MetricVars<data_t> &metric_vars,
            const vars_t<Tensor<1,data_t>> &d1,
            const diff2_vars_t<Tensor<2,data_t>> &d2
        ) const ;

};

#include "L3.impl.hpp"
#endif //L3_SIMP_H_INCLUDED