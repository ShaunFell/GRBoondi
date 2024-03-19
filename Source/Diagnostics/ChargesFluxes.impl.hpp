

#if !defined(DIAGNOSTIC_H_INCLUDED)
#error "This file should only be included through ChargesFluxes.hpp"
#endif

#ifndef DIAGNOSTIC_IMPL_H_INCLUDED
#define DIAGNOSTIC_IMPL_H_INCLUDED



template <class matter_t, class background_t>
template <class data_t>
void ChargesFluxes<matter_t, background_t>::compute(Cell<data_t> current_cell) const 
{
    /*
        See https://arxiv.org/pdf/2104.13420.pdf for conservation equations
    */

    Coordinates<data_t> coords(current_cell, m_dx, m_center);

    //load background variables
    MetricVars<data_t> metric_vars;

    //compute background
    m_background.compute_metric_background(metric_vars, coords);

    //load variables from Chombo grid and compute derivatives
    const MatterVars<data_t> matter_vars { current_cell.template load_vars<MatterVars>() };
    const MatterVars<Tensor<1,data_t>> matter_vars_d1 { m_deriv.template diff1<MatterVars>(current_cell) };
    const MatterDiff2Vars<Tensor<2,data_t>> matter_vars_d2 { m_deriv.template diff2<MatterDiff2Vars>(current_cell) };
    const MatterVars<data_t> advec { m_deriv.template advection<MatterVars>(current_cell, metric_vars.shift) };
    
    //calculate contravariant spatial metric
    const auto gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);

    //compute Christoffel symbols and spatial metric determinant
    Tensor<3, data_t>  chris_phys { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };
    const auto det_gamma = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
   
    //compute EM tensor
    const auto emtensor = m_matter.compute_emtensor(matter_vars, metric_vars, matter_vars_d1, matter_vars_d2, advec);
    
    /*
            compute conserved charges related to killing vectors in Kerr-Schild spacetime
    */

    //conserved energy
    data_t rho = -emtensor.rho * metric_vars.lapse;
    FOR1(i)
    {
        rho += metric_vars.shift[i] * emtensor.Si[i];
    }
    rho *= sqrt(det_gamma);


    //conserved angular momentum
    Tensor<1,data_t> ddphi;
    ddphi[0] = -coords.y;
    ddphi[1] = coords.x;
    ddphi[2] = 0.;

    data_t rhoJ = 0.;
    FOR1(i)
    {
        rhoJ += emtensor.Si[i] * ddphi[i];
    }
    rhoJ *= sqrt(det_gamma);



    //Compute fluxes
    Tensor<1,data_t> Ni_L;
    data_t R = coords.get_radius();
    Ni_L[0] = coords.x/R;
    Ni_L[1] = coords.y/R;
    Ni_L[2] = coords.z/R;

    //normalize using full metric
    data_t N_mod { 0.0 };
    FOR2(i,j){N_mod += gamma_UU[i][j] * Ni_L[i] * Ni_L[j]; };
    FOR1(i) { Ni_L[i] = Ni_L[i] / sqrt(N_mod); };

    //sphere area element
    data_t rho2 = simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
    data_t r2sintheta = sqrt(rho2)*R;
    Tensor<2,data_t> spherical_gamma = CoordinateTransformations::cartesian_to_spherical_LL(metric_vars.gamma,coords.x, coords.y, coords.z);
    data_t sqrt_det_Sigma = CoordinateTransformations::area_element_sphere(spherical_gamma);

    //dx/dphi to convert tensors to spherical polar
    Tensor<1,data_t> dxdphi;
    dxdphi[0] = -coords.y;
    dxdphi[1] = coords.x;
    dxdphi[2] = 0;

    //Energy flux density
    data_t Edot { 0.0 };
    FOR1(i)
    {
        Edot += metric_vars.lapse * Ni_L[i] * emtensor.rho *
                metric_vars.shift[i];

        FOR1(j)
        {
            Edot +=
                -Ni_L[i] * emtensor.Si[j] *
                (metric_vars.shift[i] * metric_vars.shift[j] +
                    metric_vars.lapse * metric_vars.lapse * gamma_UU[i][j]);
            FOR1(k)
            {
                Edot += Ni_L[i] * metric_vars.lapse * gamma_UU[i][j] *
                        metric_vars.shift[k] * emtensor.Sij[j][k];
            }
        }
    }

    // This factor of det_Sigma takes care of the surface element
    // The r2sintheta part is counted in the coordinate integration
    // so remove it here
    Edot *= sqrt_det_Sigma / r2sintheta;


    // Angular momentum flux density
    data_t Jdot = 0;
    FOR2(i, j)
    {
        Jdot +=
            -Ni_L[i] * metric_vars.shift[i] * emtensor.Si[j] * dxdphi[j];
        FOR1(k)
        {
            Jdot += metric_vars.lapse * emtensor.Sij[i][j] * dxdphi[j] *
                    gamma_UU[i][k] * Ni_L[k];
        }
    }
    Jdot *= sqrt_det_Sigma / r2sintheta;


    //assign to grid cell
    current_cell.store_vars(rho, c_rho);
    current_cell.store_vars(rhoJ, c_rhoJ);
    current_cell.store_vars(Edot, c_Edot);
    current_cell.store_vars(Jdot, c_Jdot);
    current_cell.store_vars(emtensor.rho, c_rhoE);
}


#endif //DIAGNOSTIC_IMPL_H_INCLUDED