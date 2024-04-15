#ifndef L2_SIMP_H_INCLUDED
#define L2_SIMP_H_INCLUDED


/*
This file calculates to the modifications to the equations of motion and energy momentum tensor coming from the second generalized Proca lagrangian

          -1      2          
L   =  ── ⋅ F  + α  ⋅ L 
 gp     4               2    2

                            ^^ We calculate this term
                       
*/



#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "TensorAlgebra.hpp"
#include "DefaultG.hpp"

/*
This class calculates to the modifications to the equations of motion and energy momentum tensor coming from the second generalized Proca lagrangian, 
which is just an arbitrary function.

NOTE: we assume G2 to be purely a function of the square of the Proca field, and not of the field strength tensor. The full theory is a function of both, however
we simplify things. Since we add virtual functions to the ProcaField class which can account for arbitrary modifications to the field equations and energy momentum 
tensor, the user can freely specify any modifications they want.

NOTE: This assumes the theory is L = -1/4 F^2 + alpha2 * L_2 with no other modifications due to higher lagrangians
*/


template <class G2 = DefaultG>
class L2
{
    public:
        struct params_t
        {
            double alpha2;
        };

    protected:

        G2 m_g2_function;
        params_t m_params;

        template <class data_t>
        using MatterVars = ADMProcaVars::MatterVars<data_t>;

        template <class data_t>
        using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

    public:

        L2() {}; //Default constructor for default initialization in matter class 

        L2(G2 a_G2_function, params_t a_params) : m_g2_function(a_G2_function), m_params{a_params} {};

    
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

};

#include "L2_simp.impl.hpp"
#endif //L2_SIMP_H_INCLUDED