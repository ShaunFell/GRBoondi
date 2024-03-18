#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

/* 
Example class to calculate initial data
This example sets the initial data to be a uniform magnetic field
*/

#include "BaseProcaField.hpp"
#include "ProcaField.hpp" //for proca parameters
#include "KerrSchild.hpp" //background class
#include "Cell.hpp"
#include "TensorAlgebra.hpp"
#include "Tensor.hpp"
#include "L2_simp.hpp"
#include "ADMProcaVars.hpp"


class Initial_Proca_Conditions
{

    public:

        struct params_t
        {
            double init_amplitude;
        };

    protected:
        

        template <class data_t>
        using MatterVars = typename ADMProcaVars::MatterVars<data_t>;

        template <class data_t>
        using MetricVars = typename ADMFixedBGVars::Vars<data_t>;

        using KerrParams = KerrSchild::params_t;
                
        const params_t m_params; //initial amplitude
        const ProcaField::params_t m_matter_params; //Proca mass and coupling
        const KerrParams m_Kerr_params; //black hole parameters
        
        KerrSchild m_background;
        double m_dx;

    public:
        Initial_Proca_Conditions(double a_dx, params_t a_params, ProcaField::params_t a_matter_params, KerrParams a_Kerr_params, KerrSchild a_background): m_dx{a_dx}, m_params{a_params}, m_matter_params{a_matter_params}, m_Kerr_params{a_Kerr_params}, m_background(a_background) {};

        template <class data_t>
        void compute(Cell<data_t> current_cell) const
        {
            //based off the initial conditions used in    http://arxiv.org/abs/1705.01544

            //location of cell
            Coordinates<data_t> coords(current_cell, m_dx, m_Kerr_params.center);
            
            //compute background variables
            MetricVars<data_t> metric_vars;
            m_background.compute_metric_background(metric_vars, coords);
            
            //flush all variables on cell
            MatterVars<data_t> mattervars;
            VarsTools::assign(mattervars,0.);
            
            const data_t kerrMass = m_Kerr_params.mass;
            const data_t kerrSpin = m_Kerr_params.spin;
            const data_t kerrSpin2 = kerrSpin*kerrSpin;      

            const data_t rho = coords.get_radius(); //x^2 + y^2 + z^2
            const data_t rho2 = rho * rho;

            const data_t coords_z = coords.z;

            //Calculate conformal factor
            data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
            chi = pow(chi, -1.0 / 3.0);


            // the Kerr Schild radius r
            const data_t r2 = 0.5 * (rho2 - kerrSpin2) +
                                sqrt(0.25 * (rho2 - kerrSpin2) * (rho2 - kerrSpin2) + kerrSpin2 * coords_z*coords_z);
            const data_t radius = sqrt(r2);

            //initial profile
            data_t alpha = kerrMass * m_matter_params.mass;
            data_t r0_BL { 1.0 / (m_matter_params.mass * alpha) };

            //set non-zero grid variables
            mattervars.Avec[0] = m_params.init_amplitude * pow(chi, 3.) * exp(-radius / r0_BL);

            //export to grid
            current_cell.store_vars(mattervars);
            
        }
};


#endif //INITIALCONDITIONS_HPP_