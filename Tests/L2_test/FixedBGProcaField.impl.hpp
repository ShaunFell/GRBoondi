/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FIXEDBGPROCAFIELD_HPP_)
#error "This file should only be included through FixedBGProcaField.hpp"
#endif

#ifndef FIXEDBGPROCAFIELD_IMPL_HPP_
#define FIXEDBGPROCAFIELD_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> FixedBGProcaField<potential_t>::compute_emtensor(
    const vars_t<data_t> &vars, const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1, data_t>> &d1, const Tensor<2, data_t> &gamma_UU,
    const Tensor<3, data_t> &chris_phys_ULL) const
{
    emtensor_t<data_t> out;

    // Getting all terms associated with the potential
    data_t rho_potential = 0;
    Tensor<1, data_t> Si_potential;
    Tensor<2, data_t> Sij_potential;
    m_potential.compute_stress_energy(rho_potential, Si_potential,
                                      Sij_potential, vars, d1, gamma_UU,
                                      metric_vars);

    // Some useful quantities

    // D_i A_j
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] += -chris_phys_ULL[k][i][j] * vars.Avec[k]; }
    }

    double d1Asum{0.};
    FOR2(i, j) { d1Asum += d1.Avec[i][j]; };

    double chris_sum{0.};
    FOR3(i, j, k) { chris_sum += chris_phys_ULL[i][j][k]; };

    double DASum{0.};
    FOR2(i, j) { DASum += DA[i][j]; }
    std::cout << "test DASum: " << DASum << std::endl;
    std::cout << "test d1Asum: " << d1Asum << std::endl;
    std::cout << "test chris_sum: " << chris_sum << std::endl;

    // D_i A_j - D_j A_i (NB exterior derivative, so christoffel symbols cancel)
    Tensor<2, data_t> diff_DA;
    FOR2(i, j) { diff_DA[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] = Sij_potential[i][j];

        FOR2(k, l)
        {
            out.Sij[i][j] += gamma_UU[k][l] * (diff_DA[i][k] * diff_DA[j][l]) -
                             metric_vars.gamma[i][k] * metric_vars.gamma[j][l] *
                                 vars.Evec[k] * vars.Evec[l] +
                             0.5 * metric_vars.gamma[k][l] *
                                 metric_vars.gamma[i][j] * vars.Evec[k] *
                                 vars.Evec[l];
            FOR2(m, n)
            {
                out.Sij[i][j] += -0.5 * metric_vars.gamma[i][j] *
                                 gamma_UU[k][m] * gamma_UU[l][n] * DA[k][l] *
                                 diff_DA[m][n];
            }
        }
    }

    // S = Tr_S_ij
    out.S = 0.0;
    FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; }

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        out.Si[i] = Si_potential[i];

        FOR1(j) { out.Si[i] += vars.Evec[j] * diff_DA[i][j]; }
    }

    // rho = n^a n^b T_ab
    out.rho = rho_potential;
    FOR2(i, j)
    {
        out.rho += 0.5 * metric_vars.gamma[i][j] * vars.Evec[i] * vars.Evec[j];
        FOR2(k, l)
        {
            out.rho += 0.5 * gamma_UU[i][k] * gamma_UU[j][l] * DA[k][l] *
                       diff_DA[i][j];
        }
    }

    double Enorm{0.};
    FOR2(i, j)
    {
        Enorm += vars.Evec[i] * vars.Evec[j] * metric_vars.gamma[i][j];
    }
    std::cout << "test Enorm: " << Enorm << std::endl;

    std::cout << "test rho potential: " << rho_potential << std::endl;

    return out;
}

// Adds VF evolution to the RHS
template <class potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FixedBGProcaField<potential_t>::matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const MetricVars<data_t> &metric_vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    // calculate full spatial christoffel symbols
    using namespace TensorAlgebra;
    const auto gamma_UU = compute_inverse(metric_vars.gamma);
    const auto chris_phys = compute_christoffel(metric_vars.d1_gamma, gamma_UU);

    // compute terms which are affected by potential
    // dphidt = phi time derivative excluding shift advection terms
    // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
    data_t dVdA = 0;
    data_t dphidt = 0;
    m_potential.compute_potential(dVdA, dphidt, vars, d1, metric_vars);

    // evolution equations for vector fields phi, A_i (note indices down) and
    // the conjugate momentum E^i (index up)
    total_rhs.phi = dphidt + advec.phi;

    FOR1(i)
    {
        total_rhs.Avec[i] = -metric_vars.lapse * d1.phi[i] -
                            vars.phi * metric_vars.d1_lapse[i] + advec.Avec[i];
        FOR1(j)
        {
            total_rhs.Avec[i] +=
                -metric_vars.lapse * metric_vars.gamma[i][j] * vars.Evec[j] +
                vars.Avec[j] * metric_vars.d1_shift[j][i];
        }
    }

    // variable for term (D_i A_j - D_j A_i)
    // NB Christoffel symbols cancel and take care with
    // indices - the second index is the derivative index
    Tensor<2, data_t> diff_DA;
    FOR2(i, j) { diff_DA[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; }

    // NB This is for E^i with indices up
    FOR1(i)
    {
        total_rhs.Evec[i] =
            metric_vars.lapse * metric_vars.K * vars.Evec[i] + advec.Evec[i];
        FOR1(j)
        {
            // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
            total_rhs.Evec[i] +=
                gamma_UU[i][j] * (metric_vars.lapse * d1.Z[j] +
                                  metric_vars.lapse * dVdA * vars.Avec[j]) -
                vars.Evec[j] * metric_vars.d1_shift[i][j];
        }

        FOR3(j, k, l)
        {
            total_rhs.Evec[i] +=
                gamma_UU[k][j] * gamma_UU[i][l] *
                (metric_vars.d1_lapse[k] * diff_DA[l][j] +
                 metric_vars.lapse * (d2.Avec[j][l][k] - d2.Avec[l][j][k]));

            FOR1(m)
            {
                total_rhs.Evec[i] += -gamma_UU[k][j] * gamma_UU[i][l] *
                                     metric_vars.lapse *
                                     (chris_phys.ULL[m][k][l] * diff_DA[m][j] +
                                      chris_phys.ULL[m][k][j] * diff_DA[l][m]);
            }
        }
    }

    // evolution equation for the damping term Z
    // dVdA = mu^2 ( 1 + 4 c4 A^k A_k - 12 c4 phi^2)
    // (ie the second part of the constraint, eqn 27)
    total_rhs.Z =
        metric_vars.lapse * (dVdA * vars.phi - m_vector_damping * vars.Z) +
        advec.Z;

    FOR1(i)
    {
        total_rhs.Z += metric_vars.lapse * d1.Evec[i][i];
        FOR1(j)
        {
            total_rhs.Z +=
                metric_vars.lapse * chris_phys.ULL[i][i][j] * vars.Evec[j];
        }
    }
}

#endif /* FIXEDBGPROCAFIELD_IMPL_HPP_ */
