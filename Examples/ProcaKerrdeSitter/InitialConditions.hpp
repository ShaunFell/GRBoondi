#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

/*
Example class to calculate initial data
This example sets the initial data to be a uniform magnetic field
*/

#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "Cell.hpp"
#include "KerrdeSitter.hpp" //background class
#include "L2_simp.hpp"
#include "ProcaField.hpp" //for proca parameters
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

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

    using KdSParams = KerrdeSitter::params_t;

    const params_t m_params;                    // initial amplitude
    const ProcaField::params_t m_matter_params; // Proca mass and coupling
    const KdSParams m_Kerr_params;              // black hole parameters

    KerrdeSitter m_background;
    double m_dx;

  public:
    Initial_Proca_Conditions(double a_dx, params_t a_params,
                             ProcaField::params_t a_matter_params,
                             KdSParams a_Kerr_params, KerrdeSitter a_background)
        : m_dx{a_dx}, m_params{a_params}, m_matter_params{a_matter_params},
          m_Kerr_params{a_Kerr_params}, m_background(a_background) {};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // based off the initial conditions used in
        // http://arxiv.org/abs/1705.01544

        // location of cell
        Coordinates<data_t> coords(current_cell, m_dx, m_Kerr_params.center);

        // compute background variables
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, coords);

        // flush all variables on cell
        MatterVars<data_t> mattervars;
        VarsTools::assign(mattervars, 0.);

        const data_t kerrMass = m_Kerr_params.mass;
        const data_t kerrSpin = m_Kerr_params.spin;
        const data_t kerrSpin2 = kerrSpin * kerrSpin;

        const data_t radius = coords.get_radius(); // x^2 + y^2 + z^2

        // Calculate conformal factor
        data_t gamma_det =
            TensorAlgebra::compute_determinant_sym(metric_vars.gamma);

        // initial profile
        data_t alpha = kerrMass * m_matter_params.mass;
        data_t r0{1.0 / (m_matter_params.mass * alpha)};

        // set non-zero grid variables
        mattervars.Avec[0] =
            m_params.init_amplitude * exp(-radius / r0) / gamma_det;

        // export to grid
        current_cell.store_vars(mattervars);
    }
};

#endif // INITIALCONDITIONS_HPP_