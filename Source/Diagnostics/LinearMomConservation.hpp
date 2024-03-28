/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef LINEARMOMCONSERVATION_HPP_
#define LINEARMOMCONSERVATION_HPP_

#include "ADMFixedBGVars.hpp"
#include "ADMProcaVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Calculates the momentum flux S_i with type matter_t and writes it to the
//! grid, see https://arxiv.org/pdf/2104.13420.pdf for details
template <class matter_t, class background_t> 
class LinearMomConservation
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = ADMProcaVars::MatterVars<data_t>;

    //  Need d2 of certain matter vars
    template <class data_t>
    using MatterDiff2Vars = ADMProcaVars::Diff2MatterVars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> 
    using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;                        //!< The matter object
    const double m_dx;                              //!< The grid spacing
    const background_t m_background;                //!< The metric background
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const int m_dir;

  public:
    LinearMomConservation(matter_t a_matter, background_t a_background,
                          int a_dir, double a_dx,
                          std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center), m_dir(a_dir)
    {
        if (m_dir >= GR_SPACEDIM || m_dir < 0)
        {
            MayDay::Error("The direction must be 0(x), 1(y) or 2(z)");
        }
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        

        // get the metric vars from the background
        Coordinates<data_t> coords(current_cell, m_dx, m_center);
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);
        
        // copy data from chombo gridpoint into local variables, and calculate the derivatives
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto d2 = m_deriv.template diff2<MatterDiff2Vars>(current_cell);
        const auto advec = m_deriv.template advection<MatterVars>(
            current_cell, metric_vars.shift);

        // some useful quantities
        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);

        //compute the EM tensor
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1,d2,advec);

        const auto det_gamma = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        
        //sphere area element
        const data_t rho2 = simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
        const data_t R { coords.get_radius() };
        const data_t r2sintheta = sqrt(rho2)*R;
        Tensor<2,data_t> spherical_gamma = CoordinateTransformations::cartesian_to_spherical_LL(metric_vars.gamma,coords.x, coords.y, coords.z);
        const data_t sqrt_det_Sigma = CoordinateTransformations::area_element_sphere(spherical_gamma);

        // the unit coordinate one form in the radial direction
        Tensor<1, data_t> si_L;
        si_L[0] = coords.x / R;
        si_L[1] = coords.y / R;
        si_L[2] = coords.z / R;

        //Normalize the vector using the full metric
        data_t S_mod { 0. };
        FOR2(i,j){ S_mod += gamma_UU[i][j] * si_L[i] * si_L[j]; }
        //Apply normalization
        FOR1(i) { si_L[i] *= 1. / sqrt(S_mod); }

        // see eqn (12) in https://arxiv.org/pdf/2104.13420.pdf
        data_t rhoLinMom = emtensor.Si[m_dir] * sqrt(det_gamma);

        // eqn (13)
        data_t fluxLinMom = 0.0;
        FOR1(i)
        {
            fluxLinMom += -metric_vars.shift[i] * si_L[i] * emtensor.Si[m_dir];
            FOR1(j)
            {
                fluxLinMom += metric_vars.lapse * gamma_UU[i][j] *
                              emtensor.Sij[m_dir][j] * si_L[i];
            }
        }

        // Add the volume element, taking into account that r^2 sin(theta) is 
        // already applied in the spherical surface extraction
        fluxLinMom *=  sqrt_det_Sigma / r2sintheta;

        // now the source as in eqn (19)
        data_t sourceLinMom = -emtensor.rho * metric_vars.d1_lapse[m_dir];
        FOR1(i)
        {
            sourceLinMom += emtensor.Si[i] * metric_vars.d1_shift[i][m_dir];
            FOR2(j, k)
            {
                sourceLinMom += metric_vars.lapse * gamma_UU[i][k] *
                                emtensor.Sij[k][j] *
                                chris_phys.ULL[j][i][m_dir];
            }
        }
        // add in volume factor
        sourceLinMom += sqrt(det_gamma);

        current_cell.store_vars(rhoLinMom, c_rhoLinMom);
        current_cell.store_vars(fluxLinMom, c_fluxLinMom);
        current_cell.store_vars(sourceLinMom, c_sourceLinMom);
    }
};

#endif /* LINEARMOMCONSERVATION_HPP_ */
