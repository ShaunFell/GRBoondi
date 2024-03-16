/*
implementation file for ProcaField.hpp
*/

#if !defined(PROCAFIELD_H_INCLUDED)
#error "This file should only be included through GeneralizedProcaField.hpp"
#endif

#ifndef PROCAFIELD_IMPL_H_INCLUDED
#define PROCAFIELD_IMPL_H_INCLUDED

/* //SIMD vectorization
 template <class background_t>
emtensor_t<SIMD> BaseProcaField<background_t>::compute_emtensor(
    const SIMD_vars_t &matter_vars,          // the value of the variables
    const SIMD_metric_vars_t &metric_vars, 
    const SIMD_vars_d1_t &d1, // the 1st derivatives
    const Tensor<2, SIMD> &gamma_UU,       // the inverse metric
    const Tensor<3, SIMD> &chris_phys_ULL   // conformal christoffel symbols
) const
{
    emtensor_t<SIMD> out;

    // D_i A_j  3-covariant derivative of spatial covector
    Tensor<2, SIMD> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys_ULL[k][i][j] * matter_vars.Avec[k]; };
    };

    // D_i A_j - D_j A_i
    Tensor<2, SIMD> DA_antisym;
    FOR2(i, j) { DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; };

    // Electric Field Norm
    SIMD Enorm{0};
    FOR2(i, j) { Enorm += metric_vars.gamma[i][j] * matter_vars.Evec[i] * matter_vars.Evec[j]; };

    /////Components of EM tensor

    // Eulerian Energy //Checked. Agrees with Mathematica notebook
    out.rho = 1. / 2. * Enorm;
    FOR4(i, j, k, l)
    {
        out.rho += 1. / 2. * gamma_UU[k][i] * gamma_UU[l][j] * DA[i][j] *
                   DA_antisym[k][l];
    };

    // Eulerian Momentum //Checked. Agrees with Mathematica notebook.
    FOR1(i)
    {
        out.Si[i] = 0; // zero initialize
        FOR1(j) { out.Si[i] += matter_vars.Evec[j] * DA_antisym[i][j]; };
    };

    // Eulerian Stress //Chedked. Agrees with Mathematica notebook.
    FOR2(i, j)
    {
        out.Sij[i][j] = 0; // zero initialize

        out.Sij[i][j] += 1. / 2. * metric_vars.gamma[i][j] * Enorm;

        FOR2(l, k)
        {
            out.Sij[i][j] +=
                -metric_vars.gamma[i][l] * metric_vars.gamma[j][k] * matter_vars.Evec[l] * matter_vars.Evec[k] +
                gamma_UU[k][l] * DA_antisym[i][l] * DA_antisym[j][k];

            FOR2(m, n)
            {
                out.Sij[i][j] += -metric_vars.gamma[i][j] * gamma_UU[m][l] *
                                 gamma_UU[n][k] * DA[m][n] * DA_antisym[l][k];
            };
        };
    };

    // Eulerian Stress scalar
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; };

    //add modifications
    compute_emtensor_modification(out, matter_vars, metric_vars, d1, gamma_UU, chris_phys_ULL);

    return out;
}; */

//without SIMD vectorization
/*  template <class background_t>
emtensor_t<NSIMD> BaseProcaField<background_t>::compute_emtensor(
    const NSIMD_vars_t &matter_vars,          // the value of the variables
    const NSIMD_metric_vars_t &metric_vars, 
    const NSIMD_vars_d1_t &d1, // the 1st derivatives
    const Tensor<2, NSIMD> &gamma_UU,       // the inverse metric
    const Tensor<3, NSIMD> &chris_phys_ULL   // conformal christoffel symbols
) const
{
    emtensor_t<NSIMD> out;

    // D_i A_j  3-covariant derivative of spatial covector
    Tensor<2, NSIMD> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys_ULL[k][i][j] * matter_vars.Avec[k]; };
    };

    // D_i A_j - D_j A_i
    Tensor<2, NSIMD> DA_antisym;
    FOR2(i, j) { DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; };

    // Electric Field Norm
    NSIMD Enorm{0};
    FOR2(i, j) { Enorm += metric_vars.gamma[i][j] * matter_vars.Evec[i] * matter_vars.Evec[j]; };

    /////Components of EM tensor

    // Eulerian Energy //Checked. Agrees with Mathematica notebook
    out.rho = 1. / 2. * Enorm;
    FOR4(i, j, k, l)
    {
        out.rho += 1. / 2. * gamma_UU[k][i] * gamma_UU[l][j] * DA[i][j] *
                   DA_antisym[k][l];
    };

    // Eulerian Momentum //Checked. Agrees with Mathematica notebook.
    FOR1(i)
    {
        out.Si[i] = 0; // zero initialize
        FOR1(j) { out.Si[i] += matter_vars.Evec[j] * DA_antisym[i][j]; };
    };

    // Eulerian Stress //Chedked. Agrees with Mathematica notebook.
    FOR2(i, j)
    {
        out.Sij[i][j] = 0; // zero initialize

        out.Sij[i][j] += 1. / 2. * metric_vars.gamma[i][j] * Enorm;

        FOR2(l, k)
        {
            out.Sij[i][j] +=
                -metric_vars.gamma[i][l] * metric_vars.gamma[j][k] * matter_vars.Evec[l] * matter_vars.Evec[k] +
                gamma_UU[k][l] * DA_antisym[i][l] * DA_antisym[j][k];

            FOR2(m, n)
            {
                out.Sij[i][j] += -metric_vars.gamma[i][j] * gamma_UU[m][l] *
                                 gamma_UU[n][k] * DA[m][n] * DA_antisym[l][k];
            };
        };
    };

    // Eulerian Stress scalar
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; };

    //add modifications
    compute_emtensor_modification(out, matter_vars, metric_vars, d1, gamma_UU, chris_phys_ULL);

    return out;
}; */
 


//with SIMD vectorization
/* template <class background_t>
void BaseProcaField<background_t>::matter_rhs(
    SIMD_rhs_vars_t &total_rhs,             // RHS terms for all vars
    const SIMD_vars_t &matter_vars,                // the value of the variables
    const SIMD_metric_vars_t &metric_vars,
    const SIMD_vars_d1_t &d1,       // value of 1st derivs
    const SIMD_diff2_vars_t &d2, // 2nd derivs
    const SIMD_vars_t<SIMD> &advec // value of the beta^i d_i(var) terms
) const
{

    // calculate conformal contravariant metric and conformal christoffel
    // symbols
    const Tensor<2, SIMD> gamma_UU = TensorAlgebra::compute_inverse(metric_vars.gamma);
    const Tensor<3, SIMD> chris_phys =
        TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL;
      
    // covariant derivative of spatial part of Proca field
    Tensor<2, SIMD> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
    }


    // evolution equations for spatial part of vector field (index down). 
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

    // evolution equations for Electric vector field (index up)
    FOR1(i)
    {
        total_rhs.Evec[i] = metric_vars.lapse * metric_vars.K * matter_vars.Evec[i] + advec.Evec[i];

        FOR1(j)
        {
            total_rhs.Evec[i] += - matter_vars.Evec[j] * metric_vars.d1_shift[i][j];
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

    // evolution equation for auxiliary constraint-damping scalar field Z
    total_rhs.Z = 0.;

    //Evolution equation for phi could be determined from constraint equation
    total_rhs.phi = 0;

    //add modifications
    matter_rhs_modification(total_rhs, matter_vars, metric_vars, d1, d2, advec);
} */


//wihtout SIMD vectorization
/* template <class background_t>
void BaseProcaField<background_t>::matter_rhs(
    NSIMD_rhs_vars_t &total_rhs,             // RHS terms for all vars
    const NSIMD_vars_t &matter_vars,                // the value of the variables
    const NSIMD_metric_vars_t &metric_vars,
    const NSIMD_vars_d1_t &d1,       // value of 1st derivs
    const NSIMD_diff2_vars_t &d2, // 2nd derivs
    const NSIMD_vars_t &advec // value of the beta^i d_i(var) terms
) const
{

    // calculate conformal contravariant metric and conformal christoffel
    // symbols
    const Tensor<2, NSIMD> gamma_UU = TensorAlgebra::compute_inverse(metric_vars.gamma);
    const Tensor<3, NSIMD> chris_phys =
        TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL;
      
    // covariant derivative of spatial part of Proca field
    Tensor<2, NSIMD> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
    }


    // evolution equations for spatial part of vector field (index down). 
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

    // evolution equations for Electric vector field (index up)
    FOR1(i)
    {
        total_rhs.Evec[i] = metric_vars.lapse * metric_vars.K * matter_vars.Evec[i] + advec.Evec[i];

        FOR1(j)
        {
            total_rhs.Evec[i] += - matter_vars.Evec[j] * metric_vars.d1_shift[i][j];
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

    // evolution equation for auxiliary constraint-damping scalar field Z
    total_rhs.Z = 0.;

    //Evolution equation for phi could be determined from constraint equation
    total_rhs.phi = 0;

    //add modifications
    matter_rhs_modification(total_rhs, matter_vars, metric_vars, d1, d2, advec);
} */

#endif // PROCAFIELD_IMPL_H_INCLUDED
