/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

/*
   Base class for the Proca field. Must be inherited by the Proca field class
defined by the user

NOTE:
       We use the method of 'Curiously Recurring Template Pattern' to allow
arbitrary templated modifications to the theory
       https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

       This allows us to use virtual functions that can be specified by derived
classes. e.g. see ProcaKerrBH Example
*/
#ifndef BASEPROCAFIELDTEST_H_INCLUDED
#define BASEPROCAFIELDTEST_H_INCLUDED

#include "ADMFixedBGVars.hpp" //For metric variables
#include "ADMProcaVars.hpp"   //For matter variables
#include "CCZ4RHS.hpp"
#include "DefaultBackground.hpp"      //Minkowski background as default
#include "FourthOrderDerivatives.hpp" //For calculating derivatives
#include "Tensor.hpp"                 //For performing tensorial operations
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //For user-defined variables (e.g. see EMKerrBH)
#include "VarsTools.hpp"
#include "simd.hpp" //for SIMD operations

#include "MetricVariablesInterface.hpp"
#include "ProcaField.hpp"

template <class background_t, class modification_t> class BaseProcaFieldTest
{
  protected:
    const background_t m_background;

    struct params_t
    {
    };

  public:
    // constructor, inputs are matter params
    BaseProcaFieldTest(background_t a_background)
        : m_background{a_background} {};

    template <class data_t> using Vars = ADMProcaVars::MatterVars<data_t>;

    template <class data_t>
    using Diff2Vars = ADMProcaVars::Diff2MatterVars<data_t>;

    // we set G=0, so this doenst really matter, need it just for using
    // MatterCCZ4
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL) const
    {
        return emtensor_t<data_t>{};
    }; //!< the conformal christoffel symbol

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs,             // RHS terms for all vars
        const vars_t<data_t> &matter_vars,         // the value fo the variables
        const vars_t<Tensor<1, data_t>> &d1,       // value of 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, // 2nd derivs
        const vars_t<data_t> &advec // value of the beta^i d_i(var) terms
    ) const;
};

#include "BaseProcaFieldTest.impl.hpp"
#endif // BASEPROCAFIELDTEST_H_INCLUDED
