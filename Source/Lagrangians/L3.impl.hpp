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

    data_t phi_dot { 0.  };
    compute_phi_dot(phi_dot, matter_vars, metric_vars, d1, d2);

    //The lie derivative of the scalar part of the Proca field
    data_t LNPhi { 1. / metric_vars.lapse * ( phi_dot - advec.phi) };
    

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
    data_t modified_phi_dot { 0. };
    compute_phi_dot(modified_phi_dot, matter_vars, metric_vars, d1, d2);
    total_rhs.phi +=  compute_phi_dot;

    
};



template <class G3>
template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
void L3<G3>::compute_phi_dot(
    double &phi_dot,
    const vars_t<data_t> &matter_vars,
    const MetricVars<data_t> &metric_vars,
    const vars_t<Tensor<1,data_t>> &d1,
    const diff2_vars_t<Tensor<2,data_t>> &d2
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
    
    // covariant derivative of shift vector D_i beta^j
    //first index is shift index
    Tensor<2,data_t> DB;
    FOR2(i,j)
    {
        DB[j][i] = metric_vars.d1_shift[j][i];
        FOR1(k) { 
            DB[i][j] += chris_phys[j][i][k] * metric_vars.shift[k]; 
            }
    }

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
        shift_div += gamma_UU[i][j] * DB[i][j];
    }

    // spatial derivative of christoffel symbols
    Tensor<4,data_t> dChris;
    FOR4(i,i,k,p)
    {
        dChris[k][i][j][p] = 0;
        
        FOR(l)
        {
            dChris[k][i][j][p] += -0.5 * gammaUU[i][l] * gammaUU[k][j]*metric_vars.d1_gamma[i][j][p] * (metric_vars.d1_gamma[l][j][i] + metric_vars.d1_gamma[l][i][j] - metric_vars.d1_gamma[i][j][l] ) + 0.5 * gamma_UU[k][l] * (metric_vars.d2_gamma[l][j][i][p] + metric_vars.d2_gamma[i][l][j][p] - metric_vars.d2_gamma[i][j][l][p] );
        }
    }

    //second covariant derivative of shift
    Tensor<3,data_t> DDB;
    FOR3(i,l,k)
    {
        DDB[i][l][k] = metric_vars.d2_shift[i][l][k];
        
        FOR1(j)
        {
            DDB[i][l][k] += dChris[i][l][j][k] * metric_vars.shift[j] + chris_phys[i][l][j] * metric_vars.d1_shift[j][k] + chris_phys[i][k][j] * DB[j][l] - chris_phys[j][k][l] * DB[i][j];
        }
    }

    //Covariant derivative of extrinsic curvature trace
    Tensor<1,data_t> DK;
    FOR1(i)
    {
        DK[i] = - shift_div * metric_vars.d1_lapse[i] / metric_vars.lapse / metric_vars.lapse;
        FOR1(j)
        {
            DK[i] += DDB[j][j][i] / metric_vars.lapse;
        }
    }
    
    //Second covariant derivative of lapse
    Tensor<2,data_t> DD_lapse;
    FOR2(i,j)
    {
        DD_lapse[j][i] = metric_vars.d2_lapse[j][i];

        FOR(k)
        {
            DD_lapse[j][i] -= chris_phys[k][i][j] * metric_vars.d1_lapse[k];
        }
    }

   /*  //Covariant wave equation of lapse
    Tensor<2,data_t> DDATrace;
    FOR2(i,j)
    {
        DDATrace = gamma_UU[i][j] * DD_lapse[i][j];
    } */

    //Second covariant derivative of spatial vector
    Tensor<2, Tensor<1,data_t>> DDA;
    

    //X^i X_j D_i beta^j

    //constraint algebra term
    data_t CAlg = 4 * g_func * g_prime - 6 * metric_vars.K * metric_vars.lapse * matter_vars.phi * g_prime + 2 * metric_vars.lapse * DA_div * g_prime + 2 * matter_vars.phi * shift_div * g_prime - 8 * matter_vars.phi * matter_vars.phi * g_prime * g_prime - 8 * g_func * matter_vars.phi * matter_vars.phi * g_prime2 + 4 * metric_vars.K * metric_vars.lapse * matter_vars.phi * matter_vars.phi * matter_vars.phi * g_prime2 - 4 * metric_vars.lapse * matter_vars.phi * matter_vars.phi * DA_div * g_prime2 + 4 * B_dot_X_DX * matter_vars.phi * g_prime2;

    //computed using mathematica   
    data_t phi_dot   = 8*g_func*pow(g_prime,2)*pow(matter_vars.phi,2)*metric_vars.alpha + 2*g_func*g_prime*K*matter_vars.phi*pow(metric_vars.alpha,2) - 4*pow(g_prime,2)*K*pow(matter_vars.phi,3)*pow(metric_vars.alpha,2) + 8*g_func*g_prime2*matter_vars.phi*metric_vars.alpha*matter_vars.Avec[a]*matter_vars.Evec[a] + 8*pow(g_prime,2)*matter_vars.phi*metric_vars.alpha*matter_vars.Avec[b]*matter_vars.Evec[b] + 2*g_prime*K*pow(metric_vars.alpha,2)*matter_vars.Avec[c]*matter_vars.Evec[c] - 4*g_prime2*K*pow(matter_vars.phi,2)*pow(metric_vars.alpha,2)*matter_vars.Avec[c]*matter_vars.Evec[c] + 2*g_prime*matter_vars.phi*metric_vars.alpha*matter_vars.Evec[a]*metric_vars.d1_lapse[a];
    data_t phi_dot += -2*g_prime*matter_vars.Avec[c]*matter_vars.Evec[c]*metric_vars.d1_lapse[b]*metric_vars.shift[b] - 2*g_prime*metric_vars.alpha*matter_vars.Evec[b]*metric_vars.shift[a]*DA[a][b] - 2*g_func*g_prime*matter_vars.phi*metric_vars.alpha*DB[a][a] - 2*g_prime*metric_vars.alpha*matter_vars.Avec[b]*matter_vars.Evec[b]*DB[d][d] - 2*g_prime*metric_vars.alpha*matter_vars.Avec[d]*metric_vars.shift[b]*DE[b][d] - 2*g_prime*pow(metric_vars.alpha,2)*d1.phi[a]*d1.phi[b]*gammaUU[a][b] + 4*g_prime2*pow(matter_vars.phi,2)*pow(metric_vars.alpha,2)*d1.phi[a]*d1.phi[b]*gammaUU[a][b] + 16*pow(g_prime,2)*matter_vars.phi*metric_vars.alpha*d1.phi[a]*matter_vars.Avec[b]*gammaUU[a][b];
    data_t phi_dot +=16*g_func*g_prime2*matter_vars.phi*metric_vars.alpha*d1.phi[a]*matter_vars.Avec[b]*gammaUU[a][b] - 4*g_func*g_prime*metric_vars.alpha*DA[a][b]*gammaUU[a][b] + 4*pow(g_prime,2)*pow(matter_vars.phi,2)*pow(metric_vars.alpha,2)*DA[a][b]*gammaUU[a][b] + g_func*metric_vars.alpha*DD_lapse[b][a]*gammaUU[a][b] + 2*g_prime*pow(matter_vars.phi,2)*metric_vars.alpha*DD_lapse[b][a]*gammaUU[a][b] - 2*g_prime*matter_vars.phi*metric_vars.alpha*DA[b][c]*DB[a][c]*gammaUU[a][b] - 2*g_prime*metric_vars.alpha*d1.phi[a]*matter_vars.Avec[b]*DB[c][c]*gammaUU[a][b] + 4*g_prime*K*pow(metric_vars.alpha,2)*d1.phi[a]*matter_vars.Avec[c]*gammaUU[a][c];
    data_t phi_dot += -8*g_prime2*K*pow(matter_vars.phi,2)*pow(metric_vars.alpha,2)*d1.phi[a]*matter_vars.Avec[c]*gammaUU[a][c] - 4*g_prime2*matter_vars.phi*metric_vars.alpha*d1.phi[a]*d1.phi[b]*matter_vars.Avec[c]*metric_vars.shift[b]*gammaUU[a][c] - 2*g_prime*d1.phi[a]*matter_vars.Avec[c]*metric_vars.d1_lapse[b]*metric_vars.shift[b]*gammaUU[a][c] + 2*g_prime*K*matter_vars.phi*pow(metric_vars.alpha,2)*DA[a][c]*gammaUU[a][c] - 8*pow(g_prime,2)*matter_vars.phi*matter_vars.Avec[c]*metric_vars.shift[b]*DA[b][a]*gammaUU[a][c] - 8*g_func*g_prime2*matter_vars.phi*matter_vars.Avec[c]*metric_vars.shift[b]*DA[b][a]*gammaUU[a][c];
    data_t phi_dot += - 2*g_prime*K*metric_vars.alpha*matter_vars.Avec[c]*metric_vars.shift[b]*DA[b][a]*gammaUU[a][c] - 4*pow(g_prime,2)*matter_vars.phi*metric_vars.alpha*matter_vars.Avec[c]*metric_vars.shift[b]*DA[b][a]*gammaUU[a][c] + 4*g_prime2*K*pow(matter_vars.phi,2)*metric_vars.alpha*matter_vars.Avec[c]*metric_vars.shift[b]*DA[b][a]*gammaUU[a][c] - 8*pow(g_prime,2)*matter_vars.phi*matter_vars.Avec[b]*matter_vars.Avec[c]*DB[a][b]*gammaUU[a][c] - 2*g_prime*K*metric_vars.alpha*matter_vars.Avec[b]*matter_vars.Avec[c]*DB[a][b]*gammaUU[a][c] + 4*g_prime2*K*pow(matter_vars.phi,2)*metric_vars.alpha*matter_vars.Avec[b]*matter_vars.Avec[c]*DB[a][b]*gammaUU[a][c];
    data_t phi_dot += - 2*g_prime*matter_vars.phi*metric_vars.alpha*DA[b][c]*DB[a][b]*gammaUU[a][c] + 2*g_prime*matter_vars.Avec[c]*metric_vars.shift[b]*DA[b][a]*DB[d][d]*gammaUU[a][c] + 2*g_prime*matter_vars.Avec[b]*matter_vars.Avec[c]*DB[a][b]*DB[d][d]*gammaUU[a][c] + 2*g_prime*matter_vars.Avec[d]*metric_vars.shift[b]*DA[b][c]*DB[a][c]*gammaUU[a][d] + 2*g_prime*matter_vars.Avec[d]*metric_vars.shift[b]*DA[c][a]*DB[b][c]*gammaUU[a][d] + 8*pow(g_prime,2)*pow(matter_vars.phi,2)*matter_vars.Avec[b]*metric_vars.d1_lapse[a]*gammaUU[b][a] + 8*g_func*g_prime2*pow(matter_vars.phi,2)*matter_vars.Avec[b]*metric_vars.d1_lapse[a]*gammaUU[b][a];
    data_t phi_dot += - 2*g_prime*matter_vars.phi*matter_vars.Avec[b]*metric_vars.d1_lapse[a]*DB[c][c]*gammaUU[b][a] - 2*g_prime*metric_vars.alpha*d1.phi[b]*metric_vars.shift[a]*DA[a][c]*gammaUU[b][c] + 4*g_prime2*matter_vars.phi*pow(metric_vars.alpha,2)*matter_vars.Avec[d]*matter_vars.Evec[d]*DA[b][c]*gammaUU[b][c] + 2*g_prime*metric_vars.alpha*d1.phi[a]*metric_vars.shift[a]*DA[b][c]*gammaUU[b][c] + 8*g_prime2*matter_vars.phi*pow(metric_vars.alpha,2)*d1.phi[a]*matter_vars.Avec[d]*DA[b][c]*gammaUU[a][d]*gammaUU[b][c] - 8*g_prime2*matter_vars.phi*pow(metric_vars.alpha,2)*d1.phi[b]*matter_vars.Avec[d]*DA[c][a]*gammaUU[a][d]*gammaUU[b][c];
    data_t phi_dot += 4*g_prime*metric_vars.alpha*matter_vars.Avec[d]*metric_vars.d1_lapse[b]*DA[c][a]*gammaUU[a][d]*gammaUU[b][c] - 4*g_prime2*metric_vars.alpha*matter_vars.Avec[c]*matter_vars.Avec[d]*matter_vars.Evec[c]*metric_vars.shift[f]*DA[f][b]*gammaUU[b][d] - 8*pow(g_prime,2)*metric_vars.alpha*matter_vars.Avec[c]*matter_vars.Avec[d]*DA[b][a]*gammaUU[a][c]*gammaUU[b][d] - 8*g_func*g_prime2*metric_vars.alpha*matter_vars.Avec[c]*matter_vars.Avec[d]*DA[b][a]*gammaUU[a][c]*gammaUU[b][d] + 4*g_prime2*K*matter_vars.phi*pow(metric_vars.alpha,2)*matter_vars.Avec[c]*matter_vars.Avec[d]*DA[b][a]*gammaUU[a][c]*gammaUU[b][d];
    data_t phi_dot += 4*g_prime2*metric_vars.alpha*d1.phi[e]*matter_vars.Avec[c]*matter_vars.Avec[d]*metric_vars.shift[e]*DA[b][a]*gammaUU[a][c]*gammaUU[b][d] - 2*g_prime*pow(metric_vars.alpha,2)*DA[a][c]*DA[b][d]*gammaUU[a][c]*gammaUU[b][d] + 2*g_prime*pow(metric_vars.alpha,2)*DA[b][a]*DA[d][c]*gammaUU[a][c]*gammaUU[b][d] - 4*g_prime2*metric_vars.alpha*d1.phi[b]*matter_vars.Avec[c]*matter_vars.Avec[d]*metric_vars.shift[e]*DA[e][a]*gammaUU[a][c]*gammaUU[b][d] + 4*g_prime2*matter_vars.Avec[c]*matter_vars.Avec[d]*metric_vars.shift[e]*metric_vars.shift[f]*DA[e][a]*DA[f][b]*gammaUU[a][c]*gammaUU[b][d];
    data_t phi_dot += 2*g_prime*matter_vars.phi*pow(metric_vars.alpha,2)*DK[a]*matter_vars.Avec[c]*gammaUU[c][a] + 4*g_prime*K*matter_vars.phi*metric_vars.alpha*matter_vars.Avec[c]*metric_vars.d1_lapse[a]*gammaUU[c][a] - 4*g_prime2*K*pow(matter_vars.phi,3)*metric_vars.alpha*matter_vars.Avec[c]*metric_vars.d1_lapse[a]*gammaUU[c][a] - 2*g_prime*d1.phi[b]*matter_vars.Avec[c]*metric_vars.d1_lapse[a]*metric_vars.shift[b]*gammaUU[c][a] - 2*g_prime*matter_vars.phi*matter_vars.Avec[c]*metric_vars.shift[b]*DD_lapse[b][a]*gammaUU[c][a] + 2*g_prime*metric_vars.alpha*d1.phi[b]*matter_vars.Avec[c]*DB[a][b]*gammaUU[c][a];
    data_t phi_dot += - 2*g_prime*matter_vars.phi*metric_vars.d1_lapse[b]*metric_vars.shift[a]*DA[a][c]*gammaUU[c][b] - 8*g_func*g_prime2*matter_vars.phi*matter_vars.Avec[a]*matter_vars.Avec[c]*DB[b][a]*gammaUU[c][b] - 2*g_prime*metric_vars.alpha*matter_vars.Avec[d]*metric_vars.d1_lapse[a]*DA[b][c]*gammaUU[b][c]*gammaUU[d][a] + 4*g_prime2*pow(matter_vars.phi,2)*metric_vars.alpha*matter_vars.Avec[d]*metric_vars.d1_lapse[a]*DA[b][c]*gammaUU[b][c]*gammaUU[d][a] + 2*g_prime*matter_vars.Avec[d]*matter_vars.Avec[e]*metric_vars.shift[f]*(cd(-f)(gammaUU(-a)(-c))*DB[b][c] + cd(-f)(DB[b][c])*gammaUU(-a)(-c))*gammaUU[b][e]*gammaUU[d][a];
    data_t phi_dot += 2*g_prime*matter_vars.Avec[c]*metric_vars.shift[a]*DA[a][b]*DB[d][c]*gammaUU[d][b] - 4*g_prime2*matter_vars.phi*matter_vars.Avec[c]*matter_vars.Avec[d]*metric_vars.d1_lapse[b]*metric_vars.shift[e]*DA[e][a]*gammaUU[a][c]*gammaUU[d][b] + 2*g_prime*metric_vars.shift[a]*metric_vars.shift[b]*DA[a][c]*DA[b][d]*gammaUU[d][c] + 4*g_prime2*matter_vars.Avec[b]*matter_vars.Avec[c]*matter_vars.Avec[d]*metric_vars.shift[f]*DA[f][a]*DB[e][b]*gammaUU[a][c]*gammaUU[d][e] - 4*g_prime2*matter_vars.phi*metric_vars.alpha*matter_vars.Avec[c]*matter_vars.Avec[d]*DA[e][b]*DB[a][c]*gammaUU[a][d]*gammaUU[e][b];
    data_t phi_dot += -  4*g_prime2*matter_vars.phi*metric_vars.alpha*matter_vars.Avec[d]*metric_vars.shift[b]*DA[b][a]*DA[e][c]*gammaUU[a][d]*gammaUU[e][c] + 4*g_prime2*pow(metric_vars.alpha,2)*matter_vars.Avec[d]*matter_vars.Avec[e]*DA[c][a]*DA[f][b]*gammaUU[a][d]*gammaUU[b][e]*gammaUU[f][c] - 4*g_prime2*pow(metric_vars.alpha,2)*matter_vars.Avec[d]*matter_vars.Avec[e]*DA[b][a]*DA[f][c]*gammaUU[a][d]*gammaUU[b][e]*gammaUU[f][c] - 2*g_prime*pow(metric_vars.alpha,2)*matter_vars.Avec[d]*gammaUU[a][d]*gammaUU[b][c]*DDA[a][b][c] - 2*g_prime*matter_vars.phi*metric_vars.alpha*metric_vars.shift[a]*gammaUU[b][c]*DDA[b][a][c];
    data_t phi_dot += 2*g_prime*pow(metric_vars.alpha,2)*matter_vars.Avec[d]*gammaUU[a][d]*gammaUU[b][c]*DDA[b][c][a] + 2*g_prime*matter_vars.Avec[c]*metric_vars.shift[b]*metric_vars.shift[d]*gammaUU[a][c]*DDA[d][b][a] - 2*g_prime*matter_vars.phi*metric_vars.alpha*matter_vars.Avec[b]*gammaUU[a][c]*DDB[a][c][b];






}





#endif //L3_SIMP_IMPL_H_INCLUDED
