/*
implementation file for ProcaField.hpp
*/

#if !defined(BASEPROCAFIELD_H_INCLUDED)
#error "This file should only be included through GeneralizedProcaField.hpp"
#endif

#ifndef BASEPROCAFIELD_IMPL_H_INCLUDED
#define BASEPROCAFIELD_IMPL_H_INCLUDED

template <class background_t, class modification_t>
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
emtensor_t<data_t> BaseProcaField<background_t, modification_t>::compute_emtensor(
    const vars_t<data_t> &matter_vars, //the value fo the variables
    const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, //value of 1st derivs
    const diff2_vars_t<Tensor<2, data_t>> &d2, //2nd derivs
    const vars_t<data_t> &advec //value of the beta^i d_i(var) terms
) const
{
    emtensor_t<data_t> out;

    auto gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
    auto chris_ULL { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };

    // D_i A_j  3-covariant derivative of spatial covector
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_ULL[k][i][j] * matter_vars.Avec[k]; };
    };

    // D_i A_j - D_j A_i
    Tensor<2, data_t> DA_antisym;
    FOR2(i, j) { DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; };

    // Electric Field Norm
    data_t Enorm{0};
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

    //add modifications ala CRTP
    static_cast<const modification_t*>(this)->compute_emtensor_modification(out, matter_vars, metric_vars, d1, d2, advec);

    return out;
};


template <class background_t, class modification_t>
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
void BaseProcaField<background_t, modification_t>::matter_rhs(
    rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
    const vars_t<data_t> &matter_vars, //the value fo the variables
    const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, //value of 1st derivs
    const diff2_vars_t<Tensor<2, data_t>> &d2, //2nd derivs
    const vars_t<data_t> &advec //value of the beta^i d_i(var) terms
) const
{

    // calculate conformal contravariant metric and conformal christoffel
    // symbols
    const Tensor<2, data_t> gamma_UU = TensorAlgebra::compute_inverse(metric_vars.gamma);
    const Tensor<3, data_t> chris_phys =
        TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL;
    
    // covariant derivative of spatial part of Proca field
    Tensor<2, data_t> DA;
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
    //for base Proca, covariant divergence of Proca field vanishes.
    total_rhs.phi = metric_vars.lapse * matter_vars.phi * metric_vars.K + advec.phi;
    FOR2(i,j)
    {
        total_rhs.phi += -  gamma_UU[i][j] * ( matter_vars.Avec[j] * metric_vars.d1_lapse[i] + metric_vars.lapse * DA[i][j] );
    }

    //add modifications ala CRTP
    static_cast<const modification_t*>(this)->matter_rhs_modification(total_rhs, matter_vars, metric_vars, d1, d2, advec);
};



#endif // BASEPROCAFIELD_IMPL_H_INCLUDED
