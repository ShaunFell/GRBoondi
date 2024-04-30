/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 *
 * Modified by GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

#ifndef ADMFIXEDBGVARS_HPP_
#define ADMFIXEDBGVARS_HPP_

#include "Tensor.hpp"
#include "UserVariables.hpp"
#include "VarsTools.hpp"

/// Namespace for ADM vars for fixed BG evolution
namespace ADMFixedBGVars
{
/// Vars object for ADM vars used in FixedBG evolution
template <class data_t> struct Vars
{
    // ADM vars needed in matter only rhs (ok for Proca and SF)
    Tensor<2, data_t> K_tensor;
    data_t K;

    data_t lapse;
    Tensor<1, data_t> shift;
    Tensor<2, data_t> gamma;

    Tensor<1, data_t> d1_lapse;
    Tensor<1, Tensor<1, data_t>> d1_shift;
    Tensor<2, Tensor<1, data_t>> d1_gamma;

    // Optional second derivatives of the vars
    Tensor<2, data_t> d2_lapse;
    Tensor<1, Tensor<2, data_t>> d2_shift;
    Tensor<2, Tensor<2, data_t>> d2_gamma;

    //add enum mapping to store these on the grid
    template <typename mapping_function_t>
    void enum_mapping(mapping_function_t mapping_function)
    {
        VarsTools::define_enum_mapping(mapping_function, c_K, K);
        VarsTools::define_enum_mapping(mapping_function, c_lapse, lapse);
        VarsTools::define_enum_mapping(mapping_function, GRInterval<c_shift1, c_shift3>(), shift);
    }
};

}; // namespace ADMFixedBGVars

#endif /* ADMFIXEDBGVARS_HPP_ */
