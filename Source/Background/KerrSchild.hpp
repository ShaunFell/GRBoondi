/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef KERRSCHILD_HPP_
#define KERRSCHILD_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a Kerr Schild BH
//! https://arxiv.org/pdf/gr-qc/9805023.pdf
//! https://arxiv.org/pdf/2011.07870.pdf

class KerrSchild
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        double mass = 1.0;                      //!<< The mass of the BH
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
        double spin = 0.0;                      //!< The spin param a = J / M
        bool need_2nd_derivs = false;
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    KerrSchild(params_t a_params, double a_dx) : m_params(a_params), m_dx(a_dx)
    {
        // check this spin param is sensible
        if ((m_params.spin > m_params.mass) || (m_params.spin < -m_params.mass))
        {
            MayDay::Error(
                "The dimensionless spin parameter must be in the range "
                "-1.0 < spin < 1.0");
        }
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, coords);

        // calculate and save chi
        data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        chi = pow(chi, -1.0 / 3.0);
        current_cell.store_vars(chi, c_chi);
    }

    // Kerr Schild solution
    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Coordinates<data_t> &coords) const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid including effect of spin
        // on x direction (length contraction)
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);
        const data_t cos_theta = z / r;

        // find the H and el quantities (el decomposed into space and time)
        data_t H = M * r / (r2 + a2 * cos_theta * cos_theta);
        const Tensor<1, data_t> el = {(r * x + a * y) / (r2 + a2),
                                      (r * y - a * x) / (r2 + a2), z / r};
        const data_t el_t = 1.0;

        // Calculate the gradients in el and H
        Tensor<1, data_t> dHdx;
        Tensor<1, data_t> dltdx;
        Tensor<2, data_t> dldx;
        get_KS_derivs(dHdx, dldx, dltdx, H, coords);

        // populate ADM vars
        vars.lapse = pow(1.0 + 2.0 * H * el_t * el_t, -0.5);
        FOR2(i, j)
        {
            vars.gamma[i][j] =
                TensorAlgebra::delta(i, j) + 2.0 * H * el[i] * el[j];
        }
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(vars.gamma);
        FOR1(i)
        {
            vars.shift[i] = 0;
            FOR1(j)
            {
                vars.shift[i] += gamma_UU[i][j] * 2.0 * H * el[j] * el_t;
            }
        }

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k)
        {
            vars.d1_gamma[i][j][k] =
                2.0 * (el[i] * el[j] * dHdx[k] + H * el[i] * dldx[j][k] +
                       H * el[j] * dldx[i][k]);
        }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = -pow(vars.lapse, 3.0) * el_t *
                               (el_t * dHdx[i] + 2.0 * H * dltdx[i]);
        }

        // use the fact that shift^i = lapse^2 * shift_i
        FOR2(i, j)
        {
            vars.d1_shift[i][j] =
                2.0 * el_t * dHdx[j] * pow(vars.lapse, 2.0) * el[i] +
                4.0 * el_t * H * vars.lapse * vars.d1_lapse[j] * el[i] +
                2.0 * el_t * H * pow(vars.lapse, 2.0) * dldx[i][j] +
                2.0 * dltdx[j] * H * pow(vars.lapse, 2.0) * el[i];
        }

        // calculate the extrinsic curvature, using the fact that
        // 2 * lapse * K_ij = D_i \beta_j + D_j \beta_i - dgamma_ij dt
        // and dgamma_ij dt = 0 in chosen fixed gauge
        const auto chris_phys = compute_christoffel(vars.d1_gamma, gamma_UU);
        FOR2(i, j)
        {
            vars.K_tensor[i][j] = 0.0;
            FOR1(k)
            {
                vars.K_tensor[i][j] +=
                    vars.gamma[k][j] * vars.d1_shift[k][i] +
                    vars.gamma[k][i] * vars.d1_shift[k][j] +
                    (vars.d1_gamma[k][i][j] + vars.d1_gamma[k][j][i]) *
                        vars.shift[k];
                FOR1(m)
                {
                    vars.K_tensor[i][j] += -2.0 * chris_phys.ULL[k][i][j] *
                                           vars.gamma[k][m] * vars.shift[m];
                }
            }
            vars.K_tensor[i][j] *= 0.5 / vars.lapse;
        }
        vars.K = compute_trace(gamma_UU, vars.K_tensor);

        if (m_params.need_2nd_derivs)
        {
            compute_2nd_derivatives(vars, coords);
        }

    }

  protected:
    /// Work out the gradients of the quantities H and el appearing in the Kerr
    /// Schild solution
    template <class data_t>
    void get_KS_derivs(Tensor<1, data_t> &dHdx, Tensor<2, data_t> &dldx,
                       Tensor<1, data_t> &dltdx, const data_t &H,
                       const Coordinates<data_t> &coords) const
    {
        // black hole params - spin a
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid, and useful quantities
        Tensor<1, data_t> x_vec;
        x_vec[0] = coords.x;
        x_vec[1] = coords.y;
        x_vec[2] = coords.z;
        const data_t x = coords.x;
        const data_t y = coords.y;
        const data_t z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);
        const data_t costheta = z / r;
        const data_t costheta2 = costheta * costheta;

        using namespace TensorAlgebra;
        // derivatives of r wrt actual grid coords
        Tensor<1, data_t> drhodx;
        FOR1(i) { drhodx[i] = x_vec[i] / rho; }

        //compute drdx
        Tensor<1, data_t> drdx;
        FOR1(i)
        {
            drdx[i] =
                0.5 / r *
                (rho * drhodx[i] +
                 0.5 / sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z) *
                     (drhodx[i] * rho * (rho2 - a2) +
                      delta(i, 2) * 2.0 * a2 * z));
        }
        
        Tensor<1, data_t> dcosthetadx;
        FOR1(i) { dcosthetadx[i] = -z / r2 * drdx[i] + delta(i, 2) / r; }
     
        // work out dHdx
        FOR1(i)
        {
            dHdx[i] = H * (drdx[i] / r -
                           2.0 / (r2 + a2 * costheta2) *
                               (r * drdx[i] + a2 * costheta * dcosthetadx[i]));
        }

        // note to use convention as in rest of tensors the last index is the
        // derivative index so these are d_i l_j
        FOR1(i)
        {
            // first the el_x comp
            dldx[0][i] =
                (x_vec[0] * drdx[i] + r * delta(i, 0) + a * delta(i, 1) -
                 2.0 * r * drdx[i] * (r * x_vec[0] + a * x_vec[1]) / (r2 + a2)) /
                (r2 + a2);
            // now the el_y comp
            dldx[1][i] =
                (x_vec[1] * drdx[i] + r * delta(i, 1) - a * delta(i, 0) -
                 2.0 * r * drdx[i] * (r * x_vec[1] - a * x_vec[0]) / (r2 + a2)) /
                (r2 + a2);
            // now the el_z comp
            dldx[2][i] = -x_vec[2] * drdx[i] / r2 + delta(i, 2) / r;
        }

        // then dltdi
        FOR1(i) { dltdx[i] = 0.0; }

    }

  public:
    // compute the second derivatives of the shift vector
    // These should really only be needed for higher couplings in the Proca EOM, hence we put them in their own method
    template <class data_t, template <typename> class vars_t>
    void compute_2nd_derivatives(vars_t<data_t> &vars,
                                   const Coordinates<data_t> &coords) const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid, and useful quantities
        Tensor<1, data_t> x_vec;
        x_vec[0] = coords.x;
        x_vec[1] = coords.y;
        x_vec[2] = coords.z;
        const data_t x = coords.x;
        const data_t y = coords.y;
        const data_t z = coords.z;
        const data_t rho = coords.get_radius();
        const data_t rho2 = rho * rho;

        // the Kerr Schild radius r
        const data_t r2 = 0.5 * (rho2 - a2) +
                          sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z);
        const data_t r = sqrt(r2);
        const data_t costheta = z / r;
        const data_t costheta2 = costheta * costheta;

        //we need the lapse
        compute_metric_background(vars, coords);

        using namespace TensorAlgebra;
        // derivatives of r wrt actual grid coords
        Tensor<1, data_t> drhodx;
        FOR1(i) { drhodx[i] = x_vec[i] / rho; }

        // find the H and el quantities (el decomposed into space and time)
        data_t H = M * r / (r2 + a2 * costheta * costheta);
        const Tensor<1, data_t> el = {(r * x + a * y) / (r2 + a2),
                                      (r * y - a * x) / (r2 + a2), z / r};
        const data_t el_t = 1.0;

        //compute drdx
        Tensor<1, data_t> drdx;
        FOR1(i)
        {
            drdx[i] =
                0.5 / r *
                (rho * drhodx[i] +
                 0.5 / sqrt(0.25 * (rho2 - a2) * (rho2 - a2) + a2 * z * z) *
                     (drhodx[i] * rho * (rho2 - a2) +
                      delta(i, 2) * 2.0 * a2 * z));
        }

        
        Tensor<1, data_t> dcosthetadx;
        FOR1(i) { dcosthetadx[i] = -z / r2 * drdx[i] + delta(i, 2) / r; }

        // second derivatives of rho wrt actual grid coords
        Tensor<2,data_t> drhodxx;
        FOR2(i,j)
        {
            drhodxx[i][j] = 0.0; //zero initialize 
            drhodxx[i][j] += delta(i,j) * (rho2 - x_vec[i] * x_vec[j]) - ( 1 - delta(i,j) ) * x_vec[i] * x_vec[j];
            drhodxx[i][j] *= 1. / rho / rho / rho;
        }

        Tensor<1, data_t> dHdx;
        Tensor<2, data_t> dldx;
        Tensor<1, data_t> dltdx;
        get_KS_derivs(dHdx, dldx, dltdx, H, coords);

        Tensor<2, data_t> dltdxx;
        FOR2(i,j) { dltdxx[i][j] = 0.0; }

        //compute drdxx, using Mathematica. See associated mathematica files
        //Is there a better way of computing these?
        Tensor<2,data_t> drdxx;
        drdxx[0][0] = pow(2,-0.5)*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(-(pow(a,2)*pow(drhodx[0],2)*(pow(rho,4) - (pow(a,2) + 4*pow(z,2))*(-pow(a,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)) + pow(rho,2)*(-2*pow(a,2) - 8*pow(z,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)))) + rho*(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2))*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))*drhodxx[0][0]);
        drdxx[0][1] = pow(2,-0.5)*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(-(drhodx[0]*drhodx[1]*pow(a,2)*(pow(rho,4) - (pow(a,2) + 4*pow(z,2))*(-pow(a,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)) + pow(rho,2)*(-2*pow(a,2) - 8*pow(z,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)))) + rho*(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2))*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))*drhodxx[0][1]);
        drdxx[0][2] = pow(2,-0.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(drhodx[0]*drhodx[2] - 2*rho*drhodx[0]*(-pow(a,2) + pow(rho,2))*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2)))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5) + 2*drhodx[0]*drhodx[2]*pow(rho,2)*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5) + drhodx[0]*drhodx[2]*(-pow(a,2) + pow(rho,2))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5) - rho*drhodx[0]*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1)*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))) + rho*drhodxx[0][2] + rho*(-pow(a,2) + pow(rho,2))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5)*drhodxx[0][2]);
        drdxx[1][0] = pow(2,-0.5)*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(-(drhodx[0]*drhodx[1]*pow(a,2)*(pow(rho,4) - (pow(a,2) + 4*pow(z,2))*(-pow(a,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)) + pow(rho,2)*(-2*pow(a,2) - 8*pow(z,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)))) + rho*(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2))*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))*drhodxx[0][1]);
        drdxx[1][1] = pow(2,-0.5)*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(-(pow(a,2)*pow(drhodx[1],2)*(pow(rho,4) - (pow(a,2) + 4*pow(z,2))*(-pow(a,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)) + pow(rho,2)*(-2*pow(a,2) - 8*pow(z,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5)))) + rho*(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2))*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))*drhodxx[1][1]);
        drdxx[1][2] = pow(2,-0.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(drhodx[1]*drhodx[2] - 2*rho*drhodx[1]*(-pow(a,2) + pow(rho,2))*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2)))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5) + 2*drhodx[1]*drhodx[2]*pow(rho,2)*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5) + drhodx[1]*drhodx[2]*(-pow(a,2) + pow(rho,2))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5) - rho*drhodx[1]*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1)*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))) + rho*drhodxx[1][2] + rho*(-pow(a,2) + pow(rho,2))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5)*drhodxx[1][2]);
        drdxx[2][0] = (pow(2,-0.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(2*drhodx[0]*drhodx[2] - 4*rho*drhodx[0]*(-pow(a,2) + pow(rho,2))*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2)))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5) - 2*rho*drhodx[0]*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1)*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))) + 2*rho*drhodxx[0][2] + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5)*(-2*drhodx[0]*drhodx[2]*(pow(a,2) - 3*pow(rho,2)) + 2*rho*(-pow(a,2) + pow(rho,2))*drhodxx[0][2])))/2.;
        drdxx[2][1] = (pow(2,-0.5)*pow(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5),-0.5)*(2*drhodx[1]*drhodx[2] - 4*rho*drhodx[1]*(-pow(a,2) + pow(rho,2))*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2)))*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1.5) - 2*rho*drhodx[1]*pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-1)*(2*z*pow(a,2) + rho*drhodx[2]*(-pow(a,2) + pow(rho,2) + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),0.5))) + 2*rho*drhodxx[1][2] + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5)*(-2*drhodx[1]*drhodx[2]*(pow(a,2) - 3*pow(rho,2)) + 2*rho*(-pow(a,2) + pow(rho,2))*drhodxx[1][2])))/2.;
        drdxx[2][2] = -0.25*(pow(2,-0.5)*pow(2*rho*drhodx[2] + ((8*z*pow(a,2) + 4*rho*drhodx[2]*(-pow(a,2) + pow(rho,2)))*pow(4*pow(a,2)*pow(z,2) + pow(pow(a,2) - pow(rho,2),2),-0.5))/2.,2)*pow(-pow(a,2) + pow(rho,2) + pow(4*pow(a,2)*pow(z,2) + pow(pow(a,2) - pow(rho,2),2),0.5),-1.5)) + (pow(2,-0.5)*pow(-pow(a,2) + pow(rho,2) + pow(4*pow(a,2)*pow(z,2) + pow(pow(a,2) - pow(rho,2),2),0.5),-0.5)*(2*pow(drhodx[2],2) - (pow(8*z*pow(a,2) + 4*rho*drhodx[2]*(-pow(a,2) + pow(rho,2)),2)*pow(4*pow(a,2)*pow(z,2) + pow(pow(a,2) - pow(rho,2),2),-1.5))/4. + 2*rho*drhodxx[2][2] + pow(pow(a,4) - 2*pow(a,2)*pow(rho,2) + pow(rho,4) + 4*pow(a,2)*pow(z,2),-0.5)*(4*pow(a,2) - 2*(pow(a,2) - 3*pow(rho,2))*pow(drhodx[2],2) + 2*rho*(-pow(a,2) + pow(rho,2))*drhodxx[2][2])))/2.;

        Tensor<2,data_t> dcosthetadxx;
        FOR2(i,j)
        {
            dcosthetadxx[i][j] = 0.; //zero initialize first
            dcosthetadxx[i][j] += 2 * z / r2 / r * drdx[i] * drdx[j] - z / r2 * drdxx[i][j] - delta(i,2) / r2 * drdx[j] - delta(j,2) / r2 * drdx[i];
        }

        // Now we compute the second derivatives of H
        Tensor<2,data_t> dHdxx;
        dHdxx[0][0] = H*pow(r,-1)*pow(a2*costheta2 + pow(r,2),-2)*(6*a2*costheta2*r*(a2*pow(dcosthetadx[0],2) - pow(drdx[0],2)) + 2*pow(r,3)*(-(a2*pow(dcosthetadx[0],2)) + pow(drdx[0],2)) - 2*a2*costheta*pow(r,2)*(-6*dcosthetadx[0]*drdx[0] + r*dcosthetadxx[0][0]) - 2*pow(a2,2)*pow(costheta,3)*(2*dcosthetadx[0]*drdx[0] + r*dcosthetadxx[0][0]) + (pow(a2,2)*pow(costheta,4) - pow(r,4))*drdxx[0][0]);
        dHdxx[0][1] = H*pow(r,-1)*pow(a2*costheta2 + pow(r,2),-2)*(-2*a2*dcosthetadx[0]*(costheta*drdx[1]*(a2*costheta2 - 3*pow(r,2)) + r*dcosthetadx[1]*(-3*a2*costheta2 + pow(r,2))) + 2*drdx[0]*(r*drdx[1]*(-3*a2*costheta2 + pow(r,2)) + a2*costheta*dcosthetadx[1]*(-(a2*costheta2) + 3*pow(r,2))) - 2*a2*costheta*r*(a2*costheta2 + pow(r,2))*dcosthetadxx[0][1] + (pow(a2,2)*pow(costheta,4) - pow(r,4))*drdxx[0][1]);
        dHdxx[0][2] = H*pow(r,-1)*pow(a2*costheta2 + pow(r,2),-2)*(-2*a2*dcosthetadx[0]*(costheta*drdx[2]*(a2*costheta2 - 3*pow(r,2)) + r*dcosthetadx[2]*(-3*a2*costheta2 + pow(r,2))) + 2*drdx[0]*(r*drdx[2]*(-3*a2*costheta2 + pow(r,2)) + a2*costheta*dcosthetadx[2]*(-(a2*costheta2) + 3*pow(r,2))) - 2*a2*costheta*r*(a2*costheta2 + pow(r,2))*dcosthetadxx[0][2] + (pow(a2,2)*pow(costheta,4) - pow(r,4))*drdxx[0][2]);
        dHdxx[1][0] = dHdxx[0][1];
        dHdxx[1][1] = H*pow(r,-1)*pow(a2*costheta2 + pow(r,2),-2)*(6*a2*costheta2*r*(a2*pow(dcosthetadx[1],2) - pow(drdx[1],2)) + 2*pow(r,3)*(-(a2*pow(dcosthetadx[1],2)) + pow(drdx[1],2)) - 2*a2*costheta*pow(r,2)*(-6*dcosthetadx[1]*drdx[1] + r*dcosthetadxx[1][1]) - 2*pow(a2,2)*pow(costheta,3)*(2*dcosthetadx[1]*drdx[1] + r*dcosthetadxx[1][1]) + (pow(a2,2)*pow(costheta,4) - pow(r,4))*drdxx[1][1]);
        dHdxx[1][2] = H*pow(r,-1)*pow(a2*costheta2 + pow(r,2),-2)*(-2*a2*dcosthetadx[1]*(costheta*drdx[2]*(a2*costheta2 - 3*pow(r,2)) + r*dcosthetadx[2]*(-3*a2*costheta2 + pow(r,2))) + 2*drdx[1]*(r*drdx[2]*(-3*a2*costheta2 + pow(r,2)) + a2*costheta*dcosthetadx[2]*(-(a2*costheta2) + 3*pow(r,2))) - 2*a2*costheta*r*(a2*costheta2 + pow(r,2))*dcosthetadxx[1][2] + (pow(a2,2)*pow(costheta,4) - pow(r,4))*drdxx[1][2]);
        dHdxx[2][0] = dHdxx[0][2];
        dHdxx[2][1] = dHdxx[1][2];
        dHdxx[2][2] = H*pow(r,-1)*pow(a2*costheta2 + pow(r,2),-2)*(6*a2*costheta2*r*(a2*pow(dcosthetadx[2],2) - pow(drdx[2],2)) + 2*pow(r,3)*(-(a2*pow(dcosthetadx[2],2)) + pow(drdx[2],2)) - 2*a2*costheta*pow(r,2)*(-6*dcosthetadx[2]*drdx[2] + r*dcosthetadxx[2][2]) - 2*pow(a2,2)*pow(costheta,3)*(2*dcosthetadx[2]*drdx[2] + r*dcosthetadxx[2][2]) + (pow(a2,2)*pow(costheta,4) - pow(r,4))*drdxx[2][2]);


        //ok, now the dldxx terms.
        //First index is the Li index, second two are the derivative indicies
        //again, these were computated with mathematica
        Tensor<3, data_t> dldxx;
        dldxx[0][0][0] = pow(a2 + pow(r,2),-3)*(2*drdx[0]*(pow(a2,2) + drdx[0]*(-3*a2*r*x - y*pow(a,3) + 3*a*y*pow(r,2) + x*pow(r,3)) - pow(r,4)) + (a2*x - r*(r*x + 2*a*y))*(a2 + pow(r,2)) * drdxx[0][0]);
        dldxx[0][0][1] = pow(a2 + pow(r,2),-3)*(drdx[1]*(pow(a2,2) + 2*drdx[0]*(-3*a2*r*x - y*pow(a,3) + 3*a*y*pow(r,2) + x*pow(r,3)) - pow(r,4)) + (a2 + pow(r,2))*(-2*a*r*drdx[0] + (a2*x - r*(r*x + 2*a*y))*drdxx[0][1]));
        dldxx[0][0][2] = pow(a2 + pow(r,2),-3)*(drdx[2]*(pow(a2,2) + 2*drdx[0]*(-3*a2*r*x - y*pow(a,3) + 3*a*y*pow(r,2) + x*pow(r,3)) - pow(r,4)) + (a2*x - r*(r*x + 2*a*y))*(a2 + pow(r,2))*drdxx[0][2]);
                
        dldxx[0][1][0] = dldxx[0][0][1];
        dldxx[0][1][1] = pow(a2 + pow(r,2),-3)*(2*drdx[1]*(-2*a*r*(a2 + pow(r,2)) + drdx[1]*(-3*a2*r*x - y*pow(a,3) + 3*a*y*pow(r,2) + x*pow(r,3))) + (a2*x - r*(r*x + 2*a*y))*(a2 + pow(r,2))*drdxx[1][1]);
        dldxx[0][1][2] = pow(a2 + pow(r,2),-3)*(2*drdx[2]*(-(a*r*(a2 + pow(r,2))) + drdx[1]*(-3*a2*r*x - y*pow(a,3) + 3*a*y*pow(r,2) + x*pow(r,3))) + (a2*x - r*(r*x + 2*a*y))*(a2 + pow(r,2))*drdxx[1][2]);

        dldxx[0][2][0] = dldxx[0][0][2];
        dldxx[0][2][1] = dldxx[0][1][2];
        dldxx[0][2][2] = pow(a2 + pow(r,2),-3)*((-6*a2*r*x - 2*y*pow(a,3) + 6*a*y*pow(r,2) + 2*x*pow(r,3))*pow(drdx[2],2) + (a2*x - r*(r*x + 2*a*y))*(a2 + pow(r,2))*drdxx[2][2]);

        dldxx[1][0][0] = pow(a2 + pow(r,2),-3)*(2*drdx[0]*(2*a*r*(a2 + pow(r,2)) + drdx[0]*(-3*a2*r*y + x*pow(a,3) - 3*a*x*pow(r,2) + y*pow(r,3))) + (a2 + pow(r,2))*(2*a*r*x + a2*y - y*pow(r,2))*drdxx[0][0]);
        dldxx[1][0][1] = pow(a2 + pow(r,2),-3)*(-6*a*x*drdx[0]*drdx[1]*pow(r,2) + 2*a2*r*((a - 3*y*drdx[0])*drdx[1] + a*x*drdxx[0][1]) + 2*pow(r,3)*((a + y*drdx[0])*drdx[1] + a*x*drdxx[0][1]) - pow(r,4)*(drdx[0] + y*drdxx[0][1]) + pow(a,3)*(drdx[0]*(a + 2*x*drdx[1]) + a*y*drdxx[0][1]));
        dldxx[1][0][2] = pow(a2 + pow(r,2),-3)*(2*drdx[2]*(a*r*(a2 + pow(r,2)) + drdx[0]*(-3*a2*r*y + x*pow(a,3) - 3*a*x*pow(r,2) + y*pow(r,3))) + (a2 + pow(r,2))*(2*a*r*x + a2*y - y*pow(r,2))*drdxx[0][2]);

        dldxx[1][1][0] = dldxx[1][0][1];
        dldxx[1][1][1] = pow(a2 + pow(r,2),-3)*(2*drdx[1]*(pow(a2,2) + drdx[1]*(-3*a2*r*y + x*pow(a,3) - 3*a*x*pow(r,2) + y*pow(r,3)) - pow(r,4)) + (a2 + pow(r,2))*(2*a*r*x + a2*y - y*pow(r,2))*drdxx[1][1]);
        dldxx[1][1][2] = pow(a2 + pow(r,2),-3)*(drdx[2]*(pow(a2,2) + 2*drdx[1]*(-3*a2*r*y + x*pow(a,3) - 3*a*x*pow(r,2) + y*pow(r,3)) - pow(r,4)) + (a2 + pow(r,2))*(2*a*r*x + a2*y - y*pow(r,2))*drdxx[1][2]);

        dldxx[1][2][0] = dldxx[1][0][2];
        dldxx[1][2][1] = dldxx[1][1][2];
        dldxx[1][2][2] = pow(a2 + pow(r,2),-3)*(2*(-3*a2*r*y + x*pow(a,3) - 3*a*x*pow(r,2) + y*pow(r,3))*pow(drdx[2],2) + (a2 + pow(r,2))*(2*a*r*x + a2*y - y*pow(r,2))*drdxx[2][2]);

        dldxx[2][0][0] = z*pow(r,-3)*(2*pow(drdx[0],2) - r*drdxx[0][0]);
        dldxx[2][0][1] = z*pow(r,-3)*(2*drdx[0]*drdx[1] - r*drdxx[0][1]);
        dldxx[2][0][2] = -(pow(r,-3)*(drdx[0]*(r - 2*z*drdx[2]) + r*z*drdxx[0][2]));

        dldxx[2][1][0] = dldxx[2][0][1];
        dldxx[2][1][1] = z*pow(r,-3)*(2*pow(drdx[1],2) - r*drdxx[1][1]);
        dldxx[2][1][2] = -(pow(r,-3)*(drdx[1]*(r - 2*z*drdx[2]) + r*z*drdxx[1][2]));

        dldxx[2][2][0] = dldxx[2][0][2];
        dldxx[2][2][1] = dldxx[2][1][2];
        dldxx[2][2][2] = pow(r,-3)*(2*z*pow(drdx[2],2) - r*(2*drdx[2] + z*drdxx[2][2]));


        //Now the second derivatives of the lapse, shift, and metric

        FOR2(i,j)
        {
            vars.d2_lapse[i][j] = (pow(vars.lapse,5)*(12*(el_t*dHdx[i] + 2*H*dltdx[i])*(el_t*dHdx[j] + 2*H*dltdx[j])*pow(el_t,2) - 4*(1 + 2*H*pow(el_t,2))*(2*H*dltdx[i]*dltdx[j] + pow(el_t,2)*dHdxx[i][j] + 2*el_t*(dHdx[j]*dltdx[i] + dHdx[i]*dltdx[j] + H*dltdxx[i][j]))))/4.;
        }

        FOR3(i,j,k)
        {
            vars.d2_shift[k][i][j] = 4*el_t*dHdx[j]*vars.lapse*vars.d1_lapse[i]*el[k] + 4*H*dltdx[j]*vars.lapse*vars.d1_lapse[i]*el[k] + 4*el_t*dHdx[i]*vars.lapse*vars.d1_lapse[j]*el[k] + 4*H*dltdx[i]*vars.lapse*vars.d1_lapse[j]*el[k] + 4*el_t*H*vars.d1_lapse[i]*vars.d1_lapse[j]*el[k] + 2*dHdx[j]*dltdx[i]*el[k]*pow(vars.lapse,2) + 2*dHdx[i]*dltdx[j]*el[k]*pow(vars.lapse,2) + 2*el_t*el[k]*pow(vars.lapse,2)*dHdxx[i][j] + 4*el_t*H*vars.lapse*vars.d1_lapse[j]*dldx[k][i] + 2*el_t*dHdx[j]*pow(vars.lapse,2)*dldx[k][i] + 2*H*dltdx[j]*pow(vars.lapse,2)*dldx[k][i] + 4*el_t*H*vars.lapse*vars.d1_lapse[i]*dldx[k][j] + 2*el_t*dHdx[i]*pow(vars.lapse,2)*dldx[k][j] + 2*H*dltdx[i]*pow(vars.lapse,2)*dldx[k][j] + 2*H*el[k]*pow(vars.lapse,2)*dltdxx[i][j] + 4*el_t*H*vars.lapse*el[k]*vars.d2_lapse[i][j] + 2*el_t*H*pow(vars.lapse,2)*dldxx[k][i][j];
        }

        FOR4(i,j,k,l)
        {
            vars.d2_gamma[i][j][k][l] = 2 * (dHdxx[l][k] * el[i] * el[j] + dHdx[k] * dldx[i][l] * el[j] + dHdx[k] * el[i] * dldx[j][l] +
                                                                    dHdx[l] * dldx[i][k] * el[j] + H * dldxx[i][k][l] * el[j] + H * dldx[i][k] * dldx[j][l] +
                                                                    dHdx[l] * el[i] * dldx[j][k] + H * dldx[i][l] * dldx[j][k] + H * el[i] * dldxx[j][k][l]
                                                                    );
        }
 
        //phew! done.

    }

  public:
    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    bool check_if_excised(const Coordinates<double> &coords, double buffer = 0.97) const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;

        // work out where we are on the grid
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;

        //compute inner and outer Boyer-Lindquist horizon radii
        const double r_plus = BL_outer_horizon();
        const double r_minus = BL_inner_horizon();

        // position relative to outer horizon - 1 indicates on horizon
        // less than one is within
        const double outer_horizon = sqrt((x * x + y * y) / (2.0 * M * r_plus) + z * z / r_plus / r_plus);

        // position relative to inner horizon - 1 indicates on horizon, less
        // than 1 is within
        const double inner_horizon = sqrt((x * x + y * y) / (2.0 * M * r_minus) + z * z / r_minus / r_minus);

        bool is_excised = false;
        // value less than 1 indicates we are within the horizon
        if (outer_horizon <  buffer|| inner_horizon < 1.0 / buffer)
        {
            is_excised = true;
        }
        return is_excised;
    }

    double BL_outer_horizon() const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;
        const double r_plus = M + sqrt(M * M - a2);
        
        return r_plus;
    }

    double BL_inner_horizon() const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;
        const double r_minus = M - sqrt(M * M - a2);
        
        return r_minus;
    }
};

#endif /* KERRSCHILD_HPP_ */