/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef POTENTIAL_HPP_
#define POTENTIAL_HPP_

#include "ADMFixedBGVars.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

class Potential
{
  public:
    struct params_t
    {
        double mass;
        double self_interaction;
    };

    // class params
    params_t m_params;

    // add alias for metric vars
    template <class data_t>
    using MetricVars = typename ADMFixedBGVars::template Vars<data_t>;

  public:
    //! The constructor
    Potential(params_t a_params) : m_params(a_params) {}

    //! Set the potential function for the proca field here
    template <class data_t, template <typename> class vars_t>
    void compute_stress_energy(data_t &rho_potential,
                               Tensor<1, data_t> &Si_potential,
                               Tensor<2, data_t> &Sij_potential,
                               const vars_t<data_t> &vars,
                               const vars_t<Tensor<1, data_t>> &d1,
                               const Tensor<2, data_t> &gamma_UU,
                               const MetricVars<data_t> &metric_vars) const
    {
        // defining some useful variables to ease the load

        const double msquared = pow(m_params.mass, 2.0);
        const double c4 = m_params.self_interaction;

        data_t Asquared = 0;
        const data_t phisquared = vars.phi * vars.phi;
        FOR2(i, j) { Asquared += gamma_UU[i][j] * vars.Avec[i] * vars.Avec[j]; }

        // massive terms
        rho_potential = 0.5 * msquared * phisquared + 0.5 * msquared * Asquared;
        // Self interacting terms
        rho_potential += c4 * msquared * Asquared * Asquared +
                         2. * c4 * msquared * Asquared * phisquared -
                         3. * c4 * msquared * phisquared * phisquared;

        FOR1(i)
        {
            // massive terms
            Si_potential[i] =
                msquared * vars.phi * vars.Avec[i]
                // Self interacting terms
                + 4. * c4 * msquared * vars.phi * vars.Avec[i] * Asquared -
                4. * c4 * msquared * vars.phi * vars.Avec[i] * phisquared;
        }

        FOR2(i, j)
        {
            // massive terms
            Sij_potential[i][j] = msquared * (vars.Avec[i] * vars.Avec[j] +
                                              0.5 * metric_vars.gamma[i][j] *
                                                  vars.phi * vars.phi);

            Sij_potential[i][j] +=
                -0.5 * metric_vars.gamma[i][j] * msquared * Asquared;

            // Self interacting terms
            Sij_potential[i][j] +=
                4. * c4 * msquared * vars.Avec[i] * vars.Avec[j] * Asquared -
                c4 * msquared * metric_vars.gamma[i][j] * Asquared * Asquared -
                4. * c4 * msquared * vars.Avec[i] * vars.Avec[j] * phisquared +
                2. * c4 * msquared * metric_vars.gamma[i][j] * Asquared *
                    phisquared -
                c4 * msquared * metric_vars.gamma[i][j] * phisquared *
                    phisquared;
        }
    }

    //! Set the potential function for the proca field here
    template <class data_t, template <typename> class vars_t>
    void compute_potential(data_t &dVdA, data_t &dphidt,
                           const vars_t<data_t> &vars,
                           const vars_t<Tensor<1, data_t>> &d1,
                           const MetricVars<data_t> &metric_vars) const
    {
        // calculate full spatial christoffel symbols and gamma^ij
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        // for ease of reading
        const double c4 = m_params.self_interaction;

        // Here we are defining often used terms
        // DA[i][j] = D_i A_j
        Tensor<2, data_t> DA;
        FOR2(i, j)
        {
            DA[i][j] = d1.Avec[j][i];
            FOR1(k) { DA[i][j] += -chris_phys.ULL[k][i][j] * vars.Avec[k]; }
        }

        // DAScalar = D_i A^i
        data_t DA_scalar;
        DA_scalar = 0;
        FOR2(i, j) { DA_scalar += DA[i][j] * gamma_UU[i][j]; }

        // Xsquared = X^/mu X_/mu
        data_t Xsquared;
        Xsquared = -vars.phi * vars.phi;
        FOR2(i, j) { Xsquared += gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i]; }

        // C = 1 + 4 c4 A^k A_k - 12 c4 phi^2
        data_t C = 1.0 - 12.0 * c4 * vars.phi * vars.phi;
        FOR2(i, j)
        {
            C += 4.0 * c4 * gamma_UU[i][j] * vars.Avec[j] * vars.Avec[i];
        }

        // dVdA = mu^2 ( 1 + 4 c4 (A^k A_k - phi^2))
        dVdA = pow(m_params.mass, 2.0) * (1.0 + 4.0 * c4 * Xsquared);

        // dphidt - for now the whole thing is here since it depends mainly
        // on the form of the potential - except the advection term which is in
        // the ProcaField code
        dphidt = metric_vars.lapse / C * (1.0 + 4.0 * c4 * Xsquared) *
                     (metric_vars.K * vars.phi - DA_scalar)
                 // QUESTION: Should this be lapse * Z / C  or lapse * Z??
                 - metric_vars.lapse * vars.Z / C;
        FOR1(i)
        {
            dphidt += 8.0 * c4 * vars.phi * metric_vars.lapse / C *
                      (vars.Evec[i] * vars.Avec[i]);

            FOR1(j)
            {
                dphidt +=
                    -gamma_UU[i][j] * vars.Avec[i] * metric_vars.d1_lapse[j] +
                    8.0 * c4 * vars.phi * metric_vars.lapse / C *
                        (2.0 * vars.Avec[i] * d1.phi[j] * gamma_UU[i][j]);

                FOR2(k, l)
                {
                    dphidt +=
                        -8.0 * c4 * metric_vars.lapse / C * gamma_UU[i][k] *
                            gamma_UU[j][l] * vars.Avec[i] * vars.Avec[j] *
                            DA[k][l] +
                        8.0 * c4 * vars.phi * metric_vars.lapse / C *
                            (-metric_vars.K_tensor[i][j] * vars.Avec[k] *
                             vars.Avec[l] * gamma_UU[i][k] * gamma_UU[j][l]);
                }
            }
        }
    }
};

#endif /* POTENTIAL_HPP_ */
