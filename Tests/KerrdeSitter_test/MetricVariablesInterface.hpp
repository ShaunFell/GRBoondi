/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

/*
This function provides a unifying interface between the fixed BG vars and the
CCZ4 ones
*/

#ifndef METRICVARIABLESINTERFACE_HPP_
#define METRICVARIABLESINTERFACE_HPP_

#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "BoxLoops.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

template <class data_t> using CCZ4_Vars = CCZ4RHS<>::Vars<data_t>;

template <class data_t> using FixedBGVars = ADMFixedBGVars::Vars<data_t>;

template <class data_t>
FixedBGVars<data_t> MetricVarsInterface(CCZ4_Vars<data_t> ccz4_metric_vars)
{
    FixedBGVars<data_t> fixed_bg_vars;

    fixed_bg_vars.K = ccz4_metric_vars.K;
    fixed_bg_vars.lapse = ccz4_metric_vars.lapse;
    FOR1(i) { fixed_bg_vars.shift[i] = ccz4_metric_vars.shift[i]; };

    FOR2(i, j)
    {
        fixed_bg_vars.gamma[i][j] =
            ccz4_metric_vars.h[i][j] / ccz4_metric_vars.chi;
        fixed_bg_vars.K_tensor[i][j] =
            (1. / ccz4_metric_vars.chi) *
            (ccz4_metric_vars.A[i][j] +
             1. / 3. * ccz4_metric_vars.K * ccz4_metric_vars.h[i][j]);
    }

    return fixed_bg_vars;
}

template <class data_t>
CCZ4_Vars<data_t> MetricVarsInterface(FixedBGVars<data_t> fixed_bg_vars)
{
    CCZ4_Vars<data_t> ccz4_metric_vars;
    ccz4_metric_vars.K = fixed_bg_vars.K;
    ccz4_metric_vars.lapse = fixed_bg_vars.lapse;
    FOR1(i) { ccz4_metric_vars.shift[i] = fixed_bg_vars.shift[i]; };
    auto detgamma = TensorAlgebra::compute_determinant_sym(fixed_bg_vars.gamma);
    ccz4_metric_vars.chi = pow(detgamma, -1.0 / 3.0);

    FOR2(i, j)
    {
        ccz4_metric_vars.h[i][j] =
            fixed_bg_vars.gamma[i][j] * ccz4_metric_vars.chi;
        ccz4_metric_vars.A[i][j] =
            fixed_bg_vars.K_tensor[i][j] * ccz4_metric_vars.chi -
            1.0 / 3.0 * ccz4_metric_vars.K * ccz4_metric_vars.h[i][j];
    }

    return ccz4_metric_vars;
}

#endif // METRICVARIABLESINTERFACE_HPP_