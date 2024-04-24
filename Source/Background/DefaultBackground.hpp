/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

#ifndef DEFAULTBACKGROUND_HPP_
#define DEFAULTBACKGROUND_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"

// Class which computes a Minkowski background

class DefaultBackground
{
  public:
    //! Struct for the params of the  BH
    struct params_t
    {
        std::array<double, CH_SPACEDIM> center; //!< The center of the BH
    };

  public:
    DefaultBackground(){};

    DefaultBackground(params_t a_params, double dx) {};

    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // Set spacetime background to Minkowski in Cartesian coordinates
        MetricVars<data_t> vars;
        compute_metric_background(vars);
        current_cell.store_vars(vars);
    }

    template <class data_t> void compute_metric_background(MetricVars<data_t> &vars,
                                   const Coordinates<data_t> &coords ) const
    {
        // Set spacetime background to Minkowski in Cartesian coordinates

        vars.lapse = 1.;
        vars.K = 0.;

        FOR1(i)
        {
            vars.shift[i] = 0.;
            vars.d1_lapse[i] = 0.;

            FOR1(j)
            {
                vars.d1_shift[i][j] = 0.;
                vars.K_tensor[i][j] = 0.;
                if (i == j)
                {
                    vars.gamma[i][j] = 1.;
                }
                else
                {
                    vars.gamma[i][j] = 0.;
                }
                FOR1(k) { vars.d1_gamma[i][j][k] = 0.; }
            }
        }
    }

    virtual bool check_if_excised(const Coordinates<double> &coords,
                          const double buffer = 1.0) const
    {
        return false; // Dont ever excise
    }
};

#endif /* DEFAULTBACKGROUND_HPP_ */