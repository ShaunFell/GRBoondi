/*
GRBoondi 2024
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

#ifndef KERRDE_SITTER_HPP_
#define KERRDE_SITTER_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//Class computes the Kerr-de Sitter black hole background
// https://arxiv.org/pdf/1011.0479.pdf


class KerrdeSitter
{
    public:
        struct params_t
        {
            double mass = 1.0;                                              //!<< The mass of the BH in solar masses
            double cosmo_constant = 0.0;                          //!< The cosmological constant
            double spin = 0.0;                                                //!< The spin param a = J / M
            std::array<double, CH_SPACEDIM> center;    //!< The center of the BH
        }

        template <class data_t>
        using MetricVars = ADMFixedBGVars::Vars<data_t>;

        const params_t m_params;
        const double m_dx;

        KerrdeSitter(params_t a_params, double a_dx) : m_params(a_params), m_dx(a_dx)
        {
            //Need to check whether the BH is overspun or not
                // not a trivial check, since cosmo. constant changes the spin bounds
                // But in this constructor, its only checked once on each level method

            if (m_params.cosmo_constant <= 0.0)
            {
                MayDay::Error("The cosmological constant must be positive definite");
            }
        }



        

        template <class data_t>
        void compute(Cell<data_t> current_cell) const
        {

        }


        template <class data_t>
        void compute_metric_background(MetricVars<data_t>& metric_vars, const Coordinates<data_t>& coords) const
        {

            // get background params
            data_t CC { m_params.cosmo_constant };
            data_t spin { m_params.spin };
            data_t M { m_params.mass };
            data_t spin2 { spin * spin };

            data_t x { coords[0] };
            data_t y { coords[1] };
            data_t z { coords[2] };

            // the coordinate radius, subject to floor
            data_t r { sqrt(D_TERM(x * x, + y * y, + z * z)) };
            static const double minimum_radius = 1.0e-6;
            r = simd_max(r, minimum_radius);

            //square of radius
            data_t r2 { r * r };

            // the cylindrical radius
            data_t rrho { simd_max( sqrt(x * x + y * y), minimum_radius) };
            data_t rrho2 { rho * rho };

            // set some trig functions
            data_t cos_theta { z / r };
            data_t sin_theta { rho / r };
            data_t cos_theta2 { cos_theta * cos_theta };
            data_t sin_theta2 { sin_theta * sin_theta };

            // compute the Kerr - de Sitter metric functions
            data_t delta_theta { 1 + CC / 3 * spin2 * cos_theta2 };
            data_t delta_r { r2 - 2 * M * r + spin2 - CC / 3 * r2 * ( r2 + spin2 ) };
            data_t rho2 { r2 + spin2 * cos_theta2 };
            data_t Theta { 1 + CC / 3 * spin2 };
            data_t Lambda_r { 1 - CC / 3 * r2 };
            data_t Gamma { delta_theta * ( r2 + spin2 ) * ( r2 + spin2 ) - delta_r * spin2 * sin_theta2 };
            data_t Upsilon { delta_theta * ( r2 + spin2 ) - delta_r };


            // Compute the 3+1 metric variables

            // lapse
            data_t lapse_2 { ( delta_r * delta_theta * rho2 ) / (2 * Theta * Theta * Gamma ) };
            metric_vars.lapse = sqrt(lapse_2);

            //shift
            Tensor<1,data_t> spherical_shift;
            spherical_shift[0] = 0.;
            spherical_shift[1] = 0.;
            spherical_shift[2] = - spin * Upsilon / Gamma;

            // spatial metric
            Tensor<2,data_t> spherical_gamma;
            spherical_gamma[0][0] = rho2 / delta_r;
            spherical_gamma[1][1] = rho2 / delta_theta;
            spherical_gamma[2][2] = sin_theta2 * Gamma / (Theta * Theta * rho2);
            spherical_gamma[0][1] = 0.;
            spherical_gamma[0][2] = 0.;
            spherical_gamma[1][0] = 0.;
            spherical_gamma[1][2] = 0.;
            spherical_gamma[2][0] = 0.;
            spherical_gamma[2][1] = 0.;


            // Now the 1st derivatives

            // first derivative of the metric variables
            data_t delta_r_dr { -2 / 3 * ( 3 M  + r * ( - 3 + spin2 * CC + 2 * r2 * CC) ) };
            data_t delta_theta_dtheta { - 2 / 3 * spin2 * CC * cos_theta * sin_theta };
            data_t rho2_dr { 2 * r };
            data_t rho2_dtheta { -2 * spin2 * cos_theta * sin_theta };
            data_t Lambda_r_dr { -2 / 3 * r * CC};

            //First derivatives of lapse
            Tensor<1,data_t> dlapse_dr;
            // derivative of lapse w.r.t r
            dlapse_dr[0] = - (delta_theta / (2 * metric_vars.lapse * Theta * Theta * Gamma * Gamma )) * ( - 8 * spin2 * r * sin_theta2 * delta_r * Upsilon - delta_theta * delta_r_dr * ( r2 + spin2 ) * ( r2 + spin2 ) * 2 * rho2  + 2 * delta_r * rho2_dr * Gamma );
            // derivative of lapse w.r.t theta
            dlapse_dr[1] = ( - delta_r / (2 * Theta * Theta * Gamma * Gamma * metric_vars.lapse )) * ( 8 * spin2 * ( r2 + spin2 )  * ( r2 + spin2 )  * sin_theta * cos_theta * delta_theta * delta_theta  + 2 * spin2 * delta_r * delta_theta * sin_theta * cos_theta * ( 2 * rho2 - 4 * r2 - r * spin2 ) + 2 * spin2 * sin_theta2 * delta_r * rho2 * delta_theta_dtheta + 2 * delta_theta * Gamma * rho2_dtheta );
            // derivative of lapse w.r.t. phi
            dlapse_dr[2] = 0.;

            // First derivatives of the shift
            Tensor<1,Tensor<1,data_t>> dshift_dr;
            // Zero initialize 
            FOR2(i,j) { dshift_dr[i][j] = 0;};
            // Set the only 2 non-zero components
            // derivative w.r.t. r
            dshift_dr[2][0] = (a * delta_theta) / (Gamma * Gamma) * ( - 2 * r * delta_r * ( rho2 + r2 + spin2) + (r2 + spin2) * ( rho2 * delta_r_dr + 2 * r * delta_theta * ( r2 + spin2)));
            // derivative w.r.t. theta
            dshift_dr[2][1] = (a * delta_r ) / (Gamma * Gamma) * ( - 2 * spin2 * sin_theta * cos_theta * Upsilon + (r2 + spin2 ) * delta_theta_dtheta * ( r2 + spin2 - spin2 * sin_theta2));

            // First derivatives of the spatial metric
            Tensor<2,Tensor<1,data_t>> dgamma_dr;
            // Zero initialize
            FOR3(i,j,k) { dgamma_dr[i][j][k] = 0; }
            // derivative w.r.t. r
            dgamma_dr[0][0][0] = (1. / (delta_r * delta_r)) * ( - rho2 * delta_r_dr + delta_r * rho2_dr );
            dgamma_dr[1][1][0] = rho2_dr / delta_theta;
            dgamma_dr[2][2][0] = (1. / (Theta * Theta * rho2 * rho2)) * ( (r2 + spin2) * sin_theta2 * delta_theta * (4 * r * rho2 - (r2 + spin2) * rho2_dr) + spin2 * sin_theta2 * sin_theta2 * ( delta_r * rho2_dr - rho2 * delta_r_dr));
            // derivative w.r.t. theta
            dgamma_dr[0][0][1] = rho2_dtheta / delta_r;
            dgamma_dr[1][1][1] = (1. / (delta_theta * delta_theta)) * ( delta_theta * rho2_dtheta - rho2 * delta_theta_dtheta );
            dgamma_dr[2][2][1] = (1. / (Theta * Theta * rho2 * rho2)) * ( spin2 * sin_theta2 * sin_theta * delta_r * ( -4 * rho2 * cos_theta + sin_theta * rho2_dtheta ) + (r2 + spin2) * (r2 + spin2 ) * sin_theta * ( rho2 * ( 2 * cos_theta * delta_theta + sin_theta * delta_theta_dtheta) - sin_theta * delta_theta * rho2_dtheta ));


            // Now transform to the rectangular coordinates

            //First compute the Jacobians and their derivatives
            const Tensor<2,data_t, 3> jacobian_BL_to_Cart { spher_to_cart_jac(x, y, z) };
            const Tensor<2, data_t, 3> jacobian_Cart_to_BL { spher_to_cart_jac_inv(x, y, z) };
            const Tensor<2, Tensor<1,data_t>> jacobian_BL_to_Cart_deriv { spher_to_cart_jac_deriv(x, y, z) };
            const Tensor<2, Tensor<1, data_t>> jacobian_Cart_to_BL_deriv { spher_to_cart_jac_inv_deriv(x, y, z) };

            //Now transform the non-derivative components
            FOR1(i)
            {
                metric_vars.shift[i] = 0.;

                FOR1(j)
                {
                    metric_vars.shift[i] +=   jacobian_Cart_to_BL[i][j] * spherical_shift[j];

                    metric_vars.gamma[i][j] = 0;

                    FOR1(k)
                    {
                        metric_vars.gamma[i][j] += jacobian_BL_to_Cart[i][k] * jacobian_BL_to_Cart[j][k] * spherical_gamma[k][j];
                    }
                }
            }

            // Now transform the derivatives
            FOR1(i)
            {
                metric_vars.d1_lapse[i] = 0.;

                FOR1(j)
                {
                    metric_vars.d1_lapse[i] += jacobian_BL_to_Cart[i][j] * dlapse_dr[j];

                    metric_vars.d1_shift[i][j] = 0.;

                    FOR2(k,n)
                    {
                        metric_vars.d1_shift[i][j] += dshift_dr[k][n] * jacobian_BL_to_Cart[n][j] * jacobian_BL_to_Cart[i][k] + spherical_shift[k] * jacobian_Cart_to_BL_deriv[i][k][n] * jacobian_BL_to_Cart[n][j];

                        metric_vars.d1_gamma[i][j][k] = 0.;

                        FOR1(b)
                        {
                            metric_vars.d1_gamma[i][j][k] += spherical_gamma[i][j] * ( jacobian_BL_to_Cart_deriv[n][i][k] * jacobian_BL_to_Cart[b][j] + jacobian_BL_to_Cart[a][i] * jacobian_BL_to_Cart_deriv[b][j][k]  );
                        }

                        FOR2(a,b)
                        {
                            metric_vars.d1_gamma[i][j][k] += dgamma_dr[a][b][n] * jacobian_BL_to_Cart[n][k] * jacobian_BL_to_Cart[a][i] * jacobian_BL_to_Cart[b][j];
                        }
                    }
                }
            }

            // Now calculate the extrinsic curvature and its trace, using the transformed metric quantities
            FOR2(i,j)
            {
                metric_vars.K_tensor[i][j] = 0.;

                FOR1(k)
                {
                    metric_vars.K_tensor[i][j] += (1. / (2 * metric_vars.lapse)) * ( metric_vars.shift[k] * metric_vars.d1_gamma[i][j][k] + metric_vars.gamma[k][j] * metric_vars.d1_shift[k][i] + metric_vars.gamma[i][k] * metric_vars.d1_shift[k][j]);
                }
            }

            //compute inverse spatial metric
            auto gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
            metric_vars.K = 0.;
            FOR2(i,j)
            {
                metric_vars.K += gamma_UU[i][j] * metric_vars.K_tensor[i][j];
            }

        }

    private:
        // Jacobian matrix for spherical to cartesian transformation
        template <class data_t>
        Tensor<2,data_t, 3> spher_to_cart_jac(const data_t x, const data_t y, const data_t z)
        {
            // compute common terms
            data_t r2 { simd_max( D_TERM(x * x, + y * y, + z * z) , 1e-12) }; //set floor
            data_r r { sqrt(r2) };
            data_t rho2 { simd_max( D_TERM(x * x, + y * y), 1e-12) }; //set floor
            data_t rho { sqrt(rho2) };

            // compute jacobian
            Tensor<2,data_t,3> jacobian_out;
            jacobian_out[0][0] = x / r;
            jacobian_out[0][1] = y / r;
            jacobian_out[0][2] = z / r;
            jacobian_out[1][0] = x * z / (r2 * rho);
            jacobian_out[1][1] = y * z / (r2 * rho);
            jacobian_out[1][2] = - rho / r2;
            jacobian_out[2][0] = - y / rho2;
            jacobian_out[2][1] = x / rho2;
            jacobian_out[2][2] = 0.;

            return jacobian_out;
        }

        // Derivative of the spherical to cartesian jacobian
        template <class data_t>
        Tensor<2, Tensor<1,data_t>> spher_to_cart_jac_deriv(const data_t x, const data_t y, const data_t z)
        {
            // compute common terms
            data_t r2 { simd_max( D_TERM(x * x, + y * y, + z * z) , 1e-12) }; //set floor
            data_r r { sqrt(r2) };
            data_t rho2 { simd_max( D_TERM(x * x, + y * y), 1e-12) }; //set floor
            data_t rho { sqrt(rho2) };
            data_t x2 { x * x  };
            data_t y2 { y * y  };
            data_t z2 { z * z  };

            Tensor<2, Tensor<1,data_t>> jacobian_deriv_out;
            //derivative w.r.t. x
            jacobian_deriv_out[0][0][0] = (y2 + z2 ) / ( r * r2 );
            jacobian_deriv_out[1][0][0] = z * ( -2. * x2 * x2 - x2 * y2 + y2 * y2 + y2 * z2) / (r2 * r2 * rho2 * rho);
            jacobian_deriv_out[2][0][0] = 2 * x * y / ( rho2 * rho2 );
            jacobian_deriv_out[0][1][0] = - x * y / (r2 * r);
            jacobian_deriv_out[1][1][0] = - x * y * z * (z2 + 3. * rho2) / ( r2 * r2 * rho2 * rho);
            jacobian_deriv_out[2][1][0] = ( y2 - x2 ) / (rho2 * rho2);
            jacobian_deriv_out[0][2][0] = - x * z / (r2 * r);
            jacobian_deriv_out[1][2][0] = x * (rho2 - x2)/(r2 * r2 * rho);
            jacobian_deriv_out[2][2][0] = 0.;

            //derivative w.r.t. y
            jacobian_deriv_out[0][0][1] = - x * y / (r2 * r);
            jacobian_deriv_out[1][0][1] = - x * y * z * (z2 + 3. * rho2) / ( r2 * r2 * rho2 * rho);
            jacobian_deriv_out[2][0][1] = (y2 - x2) / (rho2 * rho2);
            jacobian_deriv_out[0][1][1] = (x2 + z2) / (r2 * r);
            jacobian_deriv_out[1][1][1] = z * ( x2 * x2 - 2. * y2 * y2 + x2 * (z2 - x2)) / ( r2 * r2 * rho2 * rho);
            jacobian_deriv_out[2][1][1] = - 2 * x * y / ( rho2 * rho2 );
            jacobian_deriv_out[0][2][1] = - y * z / (r2 * r);
            jacobian_deriv_out[1][2][1] = y * (rho2 - z2) / ( r2 * r2 * rho);
            jacobian_deriv_out[2][2][1] = 0.;

            //derivative w.r.t. z
            jacobian_deriv_out[0][0][2] = - x * z / (r2 * r);
            jacobian_deriv_out[1][0][2] = x * (rho2 - z2) / (r2 * r2 * rho);
            jacobian_deriv_out[2][0][2] = 0.;
            jacobian_deriv_out[0][1][2] = - y * z / (r2 * r);
            jacobian_deriv_out[1][1][2] = y * (rho2 - z2) / (r2 * r2 * rho);
            jacobian_deriv_out[2][1][2] = 0.;
            jacobian_deriv_out[0][2][2] = rho2 / (r2 * r);
            jacobian_deriv_out[1][2][2] = 2. * z * rho / (r2 * r2);
            jacobian_deriv_out[2][2][2] = 0.;

            return jacobian_deriv_out;

        }

        // Inverse of the spherical to cartesian jacobian
        template <class data_t>
        Tensor<2, data_t, 3> spher_to_cart_jac_inv(const data_t x, const data_t y, const data_t z)
        {
            // compute common terms
            data_t r2 { simd_max( D_TERM(x * x, + y * y, + z * z) , 1e-12) }; //set floor
            data_r r { sqrt(r2) };
            data_t rho2 { simd_max( D_TERM(x * x, + y * y), 1e-12) }; //set floor
            data_t rho { sqrt(rho2) };

            //compute inverse using TensorAlgebra methods
            Tensor<2,data_t,3> inv_jacobian_out { TensorAlgebra::compute_inverse(spher_to_cart_jac(x, y, z)) };

            return inv_jacobian_out;
        }

        // Derivative of the inverse jacobian (with respect to spherical coordinates)
        Tensor< 2, Tensor<1,data_t>> spher_to_cart_jac_inv_deriv(const data_t x, const data_t y, const data_t z)
        {
            // compute common terms
            data_t r2 { simd_max( D_TERM(x * x, + y * y, + z * z) , 1e-12) }; //set floor
            data_r r { sqrt(r2) };
            data_t rho2 { simd_max( D_TERM(x * x, + y * y), 1e-12) }; //set floor
            data_t rho { sqrt(rho2) };

            Tensor< 2, Tensor<1,data_t>> inv_jacobian_deriv_out;

            //derivative w.r.t. r;
            inv_jacobian_deriv_out[0][0][0] = 0.;
            inv_jacobian_deriv_out[1][0][0] = 0.;
            inv_jacobian_deriv_out[2][0][0] = 0.;
            inv_jacobian_deriv_out[0][1][0] = x * z / (r * rho);
            inv_jacobian_deriv_out[1][1][0] = y * z / (r * rho);
            inv_jacobian_deriv_out[2][1][0] = - rho / r;
            inv_jacobian_deriv_out[0][2][0] = - y / r;
            inv_jacobian_deriv_out[1][2][0] = x / r;
            inv_jacobian_deriv_out[2][2][0] = 0.;

            //derivative w.r.t. theta
            inv_jacobian_deriv_out[0][0][1] = x * z / (r * rho);
            inv_jacobian_deriv_out[1][0][1] = y * z / (r * rho);
            inv_jacobian_deriv_out[2][0][1] = - rho / r;
            inv_jacobian_deriv_out[0][1][1] = -x;
            inv_jacobian_deriv_out[1][1][1] = -y;
            inv_jacobian_deriv_out[2][1][1] = -z;
            inv_jacobian_deriv_out[0][2][1] = - y * z / rho;
            inv_jacobian_deriv_out[1][2][1] = x * z / rho;
            inv_jacobian_deriv_out[2][2][1] = 0.;

            //derivative w.r.t. phi
            inv_jacobian_deriv_out[0][0][2] = - y / r;
            inv_jacobian_deriv_out[1][0][2] = x / r;
            inv_jacobian_deriv_out[2][0][2] = 0.;
            inv_jacobian_deriv_out[0][1][2] = -y * z / rho;
            inv_jacobian_deriv_out[1][1][2] = x * z / rho;
            inv_jacobian_deriv_out[2][1][2] = 0.;
            inv_jacobian_deriv_out[0][2][2] = -x;
            inv_jacobian_deriv_out[1][2][2] = -y;
            inv_jacobian_deriv_out[2][2][2] = 0.;

            return inv_jacobian_deriv_out;

        }




}




#endif /* KERRDE_SITTER_HPP_ */