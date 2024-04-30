/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

/*
implementation file for ProcaField.hpp
*/

#if !defined(BASEPROCAFIELDTEST_H_INCLUDED)
#error "This file should only be included through GeneralizedProcaField.hpp"
#endif

#ifndef BASEPROCAFIELDTEST_IMPL_H_INCLUDED
#define BASEPROCAFIELDTEST_IMPL_H_INCLUDED

template <class background_t, class modification_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void BaseProcaFieldTest<background_t, modification_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs,             // RHS terms for all vars
    const vars_t<data_t> &matter_vars,         // the value fo the variables
    const vars_t<Tensor<1, data_t>> &d1,       // value of 1st derivs
    const diff2_vars_t<Tensor<2, data_t>> &d2, // 2nd derivs
    const vars_t<data_t> &advec // value of the beta^i d_i(var) terms
) const
{

    // compute non-conformal metric
    CCZ4RHS<>::Vars<data_t> ccz4_vars;
    FOR2(i, j)
    {
        ccz4_vars.h[i][j] = matter_vars.h[i][j];
        ccz4_vars.A[i][j] = matter_vars.A[i][j];
    }

    ADMFixedBGVars::Vars<data_t> fixed_bg_vars;
    fixed_bg_vars.lapse = matter_vars.lapse;
    fixed_bg_vars.K = matter_vars.K;
    FOR1(i)
    {
        fixed_bg_vars.shift[i] = matter_vars.shift[i];
        fixed_bg_vars.d1_lapse[i] = d1.lapse[i];
    };
    FOR2(i, j)
    {
        fixed_bg_vars.gamma[i][j] = matter_vars.h[i][j] / matter_vars.chi;
        fixed_bg_vars.K_tensor[i][j] =
            (matter_vars.A[i][j] +
             1.0 / 3.0 * matter_vars.K * matter_vars.h[i][j]) /
            matter_vars.chi;
        fixed_bg_vars.d1_shift[i][j] = d1.shift[i][j];
    }
    FOR3(i, j, k)
    {
        fixed_bg_vars.d1_gamma[i][j][k] =
            d1.h[i][j][k] / matter_vars.chi -
            matter_vars.h[i][j] / matter_vars.chi / matter_vars.chi * d1.chi[k];
    }

    auto h_UU = TensorAlgebra::compute_inverse_sym(matter_vars.h);
    const auto non_phys_chris =
        TensorAlgebra::compute_christoffel(d1.h, h_UU).ULL;
    const Tensor<3, data_t> chris_phys = TensorAlgebra::compute_phys_chris(
        d1.chi, matter_vars.chi, matter_vars.h, h_UU, non_phys_chris);

    // calculate conformal contravariant metric and conformal christoffel
    // symbols
    const Tensor<2, data_t> gamma_UU =
        TensorAlgebra::compute_inverse(fixed_bg_vars.gamma);

    // covariant derivative of spatial part of Proca field
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
    }

    // evolution equations for spatial part of vector field (index down)
    FOR1(i)
    {
        total_rhs.Avec[i] = -fixed_bg_vars.lapse * d1.phi[i] -
                            matter_vars.phi * d1.lapse[i] + advec.Avec[i];

        FOR1(j)
        {
            total_rhs.Avec[i] += -fixed_bg_vars.lapse *
                                     fixed_bg_vars.gamma[i][j] *
                                     matter_vars.Evec[j] +
                                 matter_vars.Avec[j] * d1.shift[j][i];
        };
    };

    // evolution equations for Electric vector field (index up)
    FOR1(i)
    {
        total_rhs.Evec[i] =
            fixed_bg_vars.lapse * fixed_bg_vars.K * matter_vars.Evec[i] +
            advec.Evec[i];

        FOR1(j) { total_rhs.Evec[i] += -matter_vars.Evec[j] * d1.shift[i][j]; }

        FOR3(j, k, l)
        {
            total_rhs.Evec[i] +=
                gamma_UU[j][k] * gamma_UU[i][l] *
                (d1.lapse[j] * (d1.Avec[k][l] - d1.Avec[l][k]) +
                 fixed_bg_vars.lapse * (d2.Avec[k][l][j] - d2.Avec[l][k][j]));

            FOR1(m)
            {
                total_rhs.Evec[i] +=
                    -fixed_bg_vars.lapse * gamma_UU[j][k] * gamma_UU[i][l] *
                    (chris_phys[m][j][l] * (d1.Avec[k][m] - d1.Avec[m][k]) +
                     chris_phys[m][j][k] * (d1.Avec[m][l] - d1.Avec[l][m]));
            };
        };
    };

    // Evolution equations for phi field totally depend on the theory, so we
    // leave is up to the user to specify them for their model
    total_rhs.phi = 0.;

    // Evolution for auxiliary Z field is also left up to the user in how they
    // want to add damping terms
    total_rhs.Z = 0.;

    static_cast<const modification_t *>(this)->matter_rhs_modification(
        total_rhs, matter_vars, fixed_bg_vars, d1, d2, advec);
};

#endif // BASEPROCAFIELDTEST_IMPL_H_INCLUDED
