#if !defined(L3_SIMP_H_INCLUDED)
#error "This file should only be included through BaseProcaFieldLevel.hpp"
#endif

#ifndef L3_SIMP_IMPL_H_INCLUDED
#define L3_SIMP_IMPL_H_INCLUDED



template <class G3>
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t>
void L3<G3>::compute_emtensor_modification(
    emtensor_t<data_t> &base_emtensor,
    const vars_t<data_t> &matter_vars,
    const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1,data_t>> &d1,
    const diff2_vars_t<Tensor<2,data_t>> &d2,
    const vars_t<data_t> &advec
) const
{
    //modify EM tensor
    Tensor<2,data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };

    data_t g_func { 0 };
    data_t g_prime { 0 };
    data_t g_prime2 { 0 };
    m_g3_function.compute_function(g_func, g_prime, g_prime2, matter_vars, metric_vars, d1, d2);

    // Precalculate commonly used terms

    // covariant derivative of spatial vector
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_ULL[k][i][j] * matter_vars.Avec[k]; };
    };

    // Covariant divergence of spatial vector
    data_t DA_div { 0. };
    FOR2(i,j)
    {
        DA_div += gamma_UU[i][j] * DA[i][j];
    }

    // contravariant spatial vector
    Tensor<1, data_t> Avec_Up { 0.,0.,0. };
    FOR2(i, j) { 
        Avec_Up[i] += gamma_UU[i][j] matter_vars.Avec[j]; 
    };

    //acceleration covector
    Tensor<1, data_t> Acc { 0.,0.,0. };
    FOR1(i)
    {
        Acc[i] += metric_vars.d1_lapse[i] / metric_vars.lapse;
    }

    //inner product of spatial vector and acceleration vector
    data_t X_dot_Acc { 0. };
    FOR2(i,j)
    {
        X_dot_Acc += gamma_UU[i][j] * matter_vars.Avec[i] * Acc[j];
    }

    // double inner product of spatial vector with extrinsic curvature  K_{ij} X^{i} X^{j}
    data_t K_dot_X_X { 0. };
    FOR4(i,j,k,l)
    {
        K_dot_X_X += gamma_UU[i][k] * gamma_UU[j][l] * metric_vars.K_tensor[i][j] * matter_vars.Avec[k] * matter_vars.Avec[l];
    }

    Tensor<1,data_t> K_dot_x { 0., 0., 0. };
    FOR3(i,j,k)
    {
        K_dot_x[i] += metric_vars.K_tensor[i][j] * matter_vars.Avec[k] * gamma_UU[j][k];
    }

    // double inner product of spatial vector with covariant derivative of spatial vector  X^{i}  X^{j}  del_i X_{j}
    data_t DX_dot_X_X { 0. };
    FOR4(i,j,k,l)
    {
        DX_dot_X_X += gamma_UU[i][k] * gamma_UU[j][l] * DA[i][j] * matter_vars.Avec[k] * matter_vars.Avec[l];
    }

    // inner product of electric vector and spatial vector  E^{i} X_{i}
    data_t E_dot_X { 0. };
    FOR1(i)
    {
        E_dot_X +=  matter_vars.Evec[i] * matter_vars.Avec[i];
    }

    // Inner product of contravariant spatial vector and covariant div of phi  X^{i} del_i phi
    data_t X_dot_dphi { 0. };
    FOR2(i,j)
    {
        X_dot_dphi += gamma_UU[i][j] * matter_vars.Avec[j] * d1.phi[i];
    }

    //Inner product of contravariant spatial vector and second index of covariatn derivative of spatial vector X^{i} del_j X_{j}
    Tensor<1,data_t> X_dot_DX2 { 0.,0.,0. };
    FOR3(i,j,k)
    {
        X_dot_DX2[i] += gamma_UU[j][k] matter_vars.Avec[j] * DX[i][k];
    }

    data_t phi_dot { compute_phi_dot(matter_vars, metric_vars, d1, d2, advec) };

    //The lie derivative of the scalar part of the Proca field
    data_t LNPhi { phi_dot / metric_vars.lapse };
    

    base_emtensor.rho += 2 * m_params.alpha3 * g_prime * ( matter_vars.phi * E_dot_X + matter_vars.phi * matter_vars.phi * matter_vars.phi * metric_vars.K - matter_vars.phi * matter_vars.phi * DA_div + DX_dot_X_X + matter_vars.phi * matter_vars.phi * X_dot_Acc);

    FOR1(i)
    {
        base_emtensor.Si[i] += 2 * m_params.alpha3 * g_prime * ( matter_vars.Avec[i] * E_dot_X - matter_vars.Avec[i] * K_dot_X_X + matter_vars.Avec[i] * matter_vars.phi * matter_vars.phi * metric_vars.K - matter_vars.Avec[i] * matter_vars.phi * DA_div + matter_vars.Avec[i] * X_dot_dphi + X_dot_DX2[i] * matter_vars.phi - matter_vars.phi * matter_vars.phi * d1.phi[i] + matter_vars.phi * matter_vars.phi * K_dot_x[i] + matter_vars.phi * matter_vars.Avec[i] * X_dot_Acc);                                                                                                                                               

        FOR1(j)
        {
            base_emtensor.Sij[i][j] += 2 * m_params.alpha3 * g_prime * ( matter_vars.phi * matter_vars.Avec[i] * matter_vars.Avec[j] * metric_vars.K + metric_vars.gamma[i][j] * matter_vars.phi * E_dot_X - 2 * metric_vars.gamma[i][j] * matter_vars.phi * K_dot_X_X - matter_vars.Avec[i] * matter_vars.Avec[j] * DA_div  + metric_vars.gamma[i][j] * matter_vars.phi * matter_vars.phi * X_dot_Acc + 2 * metric_vars.gamma[i][j] * matter_vars.phi * X_dot_dphi - metric_vars.gamma[i][j] * DX_dot_X_X + matter_vars.Avec[j] * X_dot_DX2[i]  - matter_vars.phi * matter_vars.Avec[j] * d1.phi[i] + matter_vars.Avec[i] * X_dot_DX2[j] - matter_vars.phi * matter_vars.Avec[i] * d1.phi[j] + LNPhi * (metric_vars.gamma[i][j] * matter_vars.phi * matter_vars.phi - matter_vars.Avec[i] * matter_vars.Avec[j] )  + matter_vars.phi * ( matter_vars.Avec[i] * K_dot_x[j] + matter_vars.Avec[j] * K_dot_x[i] ));
        }
    }

};



template <class G3>
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
void L3<G3>::matter_rhs_modification(
    rhs_vars_t<data_t> &total_rhs,
    const vars_t<data_t> &matter_vars,
    const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1,data_t>> &d1,
    const diff2_vars_t<Tensor<2,data_t>> &d2,
    const vars_t<data_t> &advec
) const
{
    //modify RHS
    Tensor<2,data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
    auto chris_phys { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };

    data_t g_func { 0. };
    data_t g_prime { 0 };
    data_t g_prime2 { 0 };
    m_g3_function.compute_function(g_func, g_prime, g_prime2, matter_vars, metric_vars, d1, d2);

    //Define commonly used terms and for ease of reading
    
    //del^i lapse
    Tensor<1,data_t> d_alpha_UP { 0.,0.,0. };
    FOR2(i,j)
    {
        d_alpha_UP[i] += gamma_UU[i][j] * metric_vars.d1_lapse[j];
    }

    //X^i
    Tensor<1, data_t> Avec_UP { 0.,0.,0. };
    FOR2(i,j)
    {
        Avec_UP[i] += gamma_UU[i][j] * matter_vars.Avec[j];
    }

    // covariant derivative of spatial part of Proca field
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
    }

    //Covariant divergence of spatial vector D_i X^i
    data_t DA_div { 0. };
    FOR2(i, j) { DA_div += gamma_UU[i][j] DA[i][j]; };

    // X^a del^i X_a
    Tensor<1,data_t> X_dot_DX_UP { 0.,0.,0. };
    FOR2(i,j)
    {
        X_dot_DX_UP[i] += Avec_UP[j] * gamma_UU [i][j] * DA[i][j];
    }

    //del^i phi
    Tensor<1,data_t> d_phi_UP { 0.,0.,0. };
    FOR2(i,j)
    {
        d_phi_UP[i] += gamma_UU[i][j] * metric_vars.d1_phi[j];
    }

    // Spatial vector rhs remains unchanged

    //Modify electric part
    FOR2(i,j)
    {
        total_rhs.Evec[i] += g_func * d_alpha_UP[i] - 4 * g_func * Avec_UP[i] * g_prime + 2 * metric_vars.K * Avec_UP[i] * metric_vars.lapse * matter_vars.phi * g_prime - 2 * Avec_UP[i] * metric_vars.lapse * DA_div * g_prime + 2 * Avec_UP[i] * advec.phi * g_prime + 2 * matter_vars.lapse * X_dot_DX_UP[i] * g_prime - 2 * metric_vars.lapse * matter_vars.phi * d_phi_UP[i] * g_prime;
    }

    //Spatial part remains unchanged

    //Modify scalar part

    total_rhs.phi +=  0;

    
};



template <class G3>
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
void L3<G3>::compute_phi_dot(
    const vars_t<data_t> &matter_vars,
    const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1,data_t>> &d1,
    const diff2_vars_t<Tensor<2,data_t>> &d2,
    const vars_t<data_t> &advec
) const
{
    // compute time derivative of scalar part


    Tensor<2,data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
    auto chris_phys { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };

    data_t g_func { 0. };
    data_t g_prime { 0 };
    data_t g_prime2 { 0 };
    m_g3_function.compute_function(g_func, g_prime, g_prime2, matter_vars, metric_vars, d1, d2);


    // covariant derivative of spatial part of Proca field
    Tensor<2, data_t> DA;
    FOR2(i, j)
    {
        DA[i][j] = d1.Avec[j][i];
        FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
    };
    
    // Covariant divergence of spatial vector
    data_t DA_div { 0. };
    FOR2(i,j)
    {
        DA_div += gamma_UU[i][j] * DA[i][j];
    }

    // Covariant divergence of shift vector
    data_t shift_div { 0. };
    FOR2(i,j)
    {
        shift_div += gamma_UU[i][j] * metric_vars.d1_shift[j][i] + chris_phy[i][i][j] * metric_vars.shift[j];
    }

    data_t X_dot_B_DX { 0. };
    FOR4(i,j,k,l)
    {
        B_dot_X_DX += gamma_UU[i][k] * gamma_UU[j][l] * matter_vars.Avec[i] * metric_vars.shift[j] * DA[l][k];
    }

    //constraint algebra term
    data_t CAlg = 4 * g_func * g_prime - 6 * metric_vars.K * metric_vars.lapse * matter_vars.phi * g_prime + 2 * metric_vars.lapse * DA_div * g_prime + 2 * matter_vars.phi * shift_div * g_prime - 8 * matter_vars.phi * matter_vars.phi * g_prime * g_prime - 8 * g_func * matter_vars.phi * matter_vars.phi * g_prime2 + 4 * metric_vars.K * metric_vars.lapse * matter_vars.phi * matter_vars.phi * matter_vars.phi * g_prime2 - 4 * metric_vars.lapse * matter_vars.phi * matter_vars.phi * DA_div * g_prime2 + 4 * B_dot_X_DX * matter_vars.phi * g_prime2;

}





#endif //L3_SIMP_IMPL_H_INCLUDED
