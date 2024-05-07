/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

#if !defined(L2_SIMP_H_INCLUDED)
#error "This file should only be included through BaseProcaFieldLevel.hpp"
#endif

#ifndef L2_SIMP_IMPL_H_INCLUDED
#define L2_SIMP_IMPL_H_INCLUDED

// DEBUGGING
#include "DefaultG.hpp"
#include "ProcaField.hpp"

template <class G2>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void L2<G2>::compute_emtensor_modification(
    emtensor_t<data_t> &base_emtensor, const vars_t<data_t> &matter_vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // modify EM tensor
    Tensor<2, data_t> gamma_UU{
        TensorAlgebra::compute_inverse_sym(metric_vars.gamma)};

    data_t g_func{0};
    data_t g_prime{0};
    data_t g_prime2{0};
    m_g2_function.compute_function(g_func, g_prime, g_prime2, matter_vars,
                                   metric_vars);

    // compute addition to Eulerian energy
    base_emtensor.rho +=
        m_params.alpha2 *
        (g_func + 2 * matter_vars.phi * matter_vars.phi * g_prime);

    // compute addition to Eulerian momentum and stress
    FOR1(i)
    {
        // addition to Eulerian momentum
        base_emtensor.Si[i] += m_params.alpha2 * (2 * matter_vars.Avec[i] *
                                                  matter_vars.phi * g_prime);

        FOR1(j)
        {
            // trace of addition to stress tensor
            base_emtensor.S +=
                m_params.alpha2 *
                (-3 * g_func + 2 * gamma_UU[i][j] * matter_vars.Avec[i] *
                                   matter_vars.Avec[j] * g_prime);

            // addition to stress tensor
            base_emtensor.Sij[i][j] +=
                m_params.alpha2 *
                (-metric_vars.gamma[i][j] * g_func +
                 2 * matter_vars.Avec[i] * matter_vars.Avec[j] * g_prime);
        }
    }
};

template <class G2>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void L2<G2>::matter_rhs_modification(rhs_vars_t<data_t> &total_rhs,
                                     const vars_t<data_t> &matter_vars,
                                     const MetricVars<data_t> &metric_vars,
                                     const vars_t<Tensor<1, data_t>> &d1,
                                     const diff2_vars_t<Tensor<2, data_t>> &d2,
                                     const vars_t<data_t> &advec) const
{
    // modify RHS
    Tensor<2, data_t> gamma_UU{
        TensorAlgebra::compute_inverse_sym(metric_vars.gamma)};
    auto chris_phys{
        TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL};

    // Compute G2 function
    data_t V{0.};
    data_t dVdA{0.};
    data_t dVddA{0.};
    m_g2_function.compute_function(V, dVdA, dVddA, matter_vars, metric_vars);

    // Compute constraint algebra term
    data_t gnn{dVdA - 2.0 * dVddA * matter_vars.phi * matter_vars.phi};

    // covariant derivative of spatial part of Proca field
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
    }

    // Spatial part remains unchanged

    // Modify scalar part
    // Also modify electric part in same loop
    total_rhs.phi +=
        metric_vars.lapse * dVdA * matter_vars.phi * metric_vars.K / (gnn) +
        advec.phi;
    FOR1(i)
    {
        total_rhs.phi += 2 * metric_vars.lapse * dVddA * matter_vars.phi *
                         matter_vars.Avec[i] * matter_vars.Evec[i] / gnn;

        FOR1(j)
        {
            total_rhs.phi +=
                gamma_UU[i][j] *
                (-metric_vars.lapse * dVdA / gnn * DA[i][j] -
                 matter_vars.Avec[i] * metric_vars.d1_lapse[j] +
                 2 * metric_vars.lapse * dVddA / gnn * 2 * matter_vars.phi *
                     matter_vars.Avec[i] * d1.phi[j]);

            // Modify electric part
            total_rhs.Evec[i] += 2 * metric_vars.lapse * dVdA * gamma_UU[i][j] *
                                 matter_vars.Avec[j];

            FOR2(k, l)
            {
                total_rhs.phi -=
                    gamma_UU[i][k] * gamma_UU[j][l] *
                    (2 * metric_vars.lapse * dVddA / gnn * matter_vars.phi *
                         matter_vars.Avec[i] * matter_vars.Avec[j] *
                         metric_vars.K_tensor[k][l] +
                     2 * metric_vars.lapse * dVddA / gnn * matter_vars.Avec[i] *
                         matter_vars.Avec[j] * DA[k][l]);
            }
        }
    }
};

#endif // L2_SIMP_IMPL_H_INCLUDED