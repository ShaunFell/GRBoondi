#ifndef SECONDCLASSCONSTRAINT_HPP_INCLUDED
#define SECONDCLASSCONSTRAINT_HPP_INCLUDED

#include "ProcaField.hpp"
#include "ADMProcaVars.hpp"
#include "Cell.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"


template <class background_t>
class SecondClassConstraint
{
    protected:

        const double m_dx;
        const FourthOrderDerivatives m_deriv;
        const std::array<double, CH_SPACEDIM> m_center;
        const background_t m_background;
        const ProcaField m_proca_field;


        //typename the matter variables
        template <class data_t>
        using MatterVars = typename ADMProcaVars::MatterVars<data_t>;
        
        //typename the metric variables
        template <class data_t>
        using MetricVars = typename ADMFixedBGVars::Vars<data_t>;

    public:
        SecondClassConstraint(double dx, const std::array<double, CH_SPACEDIM> a_center, const background_t a_background, const ProcaField a_proca_field):
            m_dx{dx}, m_deriv(dx), m_background(a_background), m_center(a_center), m_proca_field(a_proca_field)
        {
        };

        template <class data_t>
        void compute(Cell<data_t> current_cell) const
        {
            MetricVars<data_t> metric_vars;

            //where are we on the grid
            const Coordinates<data_t> coords(current_cell, m_dx, m_center);

            //compute derivatives
            const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);

            //compute the metric background
            m_background.compute_metric_background(metric_vars, coords);

            //load matter variables from Chombo grid
            const MatterVars<data_t> matter_vars = current_cell.template load_vars<MatterVars>();

            //compute contravariant spatial metric
            const Tensor<2,data_t> gamma_UU { TensorAlgebra::compute_inverse_sym(metric_vars.gamma) };
            auto chris_phys { TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL };

            //compute the G2 function
            data_t g_func { 0.0 };
            data_t g_prime { 0.0 };
            data_t g_prime2 { 0.0 };
            m_proca_field.m_G3.compute_function(g_func, g_prime, g_prime2, matter_vars, metric_vars);

            Tensor<2, data_t> DA;
            FOR2(i, j)
            {
                DA[i][j] = d1.Avec[j][i];
                FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
            };

            Tensor<2,data_t> DB;
            FOR2(i,j)
            {
                DB[i][j] = metric_vars.d1_shift[j][i];
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

            
            data_t B_dot_X_DX { 0. };
            FOR3(i,j,a)
            {
                B_dot_X_DX += gamma_UU[i][j] * matter_vars.Avec[i] * metric_vars.shift[a] * DA[a][j];
            }

            //compute the constraint algebra term
            data_t CAlg = 4 * g_func * g_prime - 6 * metric_vars.K * metric_vars.lapse * matter_vars.phi * g_prime + 2 * metric_vars.lapse * DA_div * g_prime + 2 * matter_vars.phi * shift_div * g_prime - 8 * matter_vars.phi * matter_vars.phi * g_prime * g_prime - 8 * g_func * matter_vars.phi * matter_vars.phi * g_prime2 + 4 * metric_vars.K * metric_vars.lapse * matter_vars.phi * matter_vars.phi * matter_vars.phi * g_prime2 - 4 * metric_vars.lapse * matter_vars.phi * matter_vars.phi * DA_div * g_prime2 + 4 * B_dot_X_DX * matter_vars.phi * g_prime2;

            //save back to chombo grid
            current_cell.store_vars(CAlg, c_gnn);
        }
};




#endif //SECONDCLASSCONSTRAINT_HPP_INCLUDED