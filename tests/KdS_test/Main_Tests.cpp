
// background includes
#include "ADMFixedBGVars.hpp"
#include "KerrdeSitter.hpp"

// GRChombo includes
#include "Cell.hpp"
#include "Coordinates.hpp"

// Chombo includes
#include "IntVect.H"
#include "UsingNamespace.H"

// c++ includes
#include <algorithm>
#include <iostream>
#include <string>

#define ERR 1e-11

int main(int argc, char *argv[])
{
    std::cout << std::string(10, '#') << std::endl;
    std::cout << "Running GRBoondi tests" << std::endl;
    std::cout << "current tests: KerrdeSitter" << std::endl;
    std::cout << std::string(10, '#') << std::endl;

    std::string tab{"\t"};

    // Here we test the outputs of the background functions
    {

        // KerrdeSitter
        {
            std::cout << "\n\n KerrdeSitter test:" << std::endl;

            double dx{1e-3};

            // initialize an IntVect object
            IntVect intvect;
            intvect[0] = 0.5;
            intvect[1] = 0.5;
            intvect[2] = 0.5;

            // create a coordinate object
            // We set the coordinates to some arbitrary value, but which
            // correspond to the
            //       values in the Mathematica notebook
            Coordinates<double> coords(
                intvect, dx,
                {0, 0, 0}); // should correspond to the point (x,y,z) = (1,1,1)
            const double x = coords.x;
            const double y = coords.y;
            const double z = coords.z;
            std::cout << tab << "(x,y,z) = (" << x << "," << y << "," << z
                      << ")" << std::endl;

            // initialize the kerr params and the kerr object
            KerrdeSitter::params_t kerr_params;
            kerr_params.mass = 1.4;
            kerr_params.spin = 0.9;
            kerr_params.cosmo_constant = 0.001;
            // kerr_params.center = {0.0, 0.0, 0.0};
            KerrdeSitter kerr_init(kerr_params, dx);
            std::cout << tab << " kerr mass = " << kerr_params.mass
                      << std::endl;
            std::cout << tab << " kerr spin = " << kerr_params.spin
                      << std::endl;
            std::cout << tab
                      << " kerr cosmo_constant = " << kerr_params.cosmo_constant
                      << std::endl;

            // Discriminant test
            double Q{KerrdeSitter::discriminant<double>(kerr_params)};
            std::cout << tab << "Q = " << Q << std::endl;
            double mathematica_result_Q{0.006037383095635};
            if (abs(Q - mathematica_result_Q) > ERR)
            {
                std::cout << tab << "Discriminant test failed" << std::endl;
            }
            else
            {
                std::cout << tab << "Discriminant test passed" << std::endl;
            }

            // Metric variable tests
            std::cout << "Metric variable tests: " << std::endl;
            // instantiate the metric variables
            ADMFixedBGVars::Vars<double> metric_vars;

            IntVect mv_intvec;
            mv_intvec[0] = 10000;
            mv_intvec[1] = 10000;
            mv_intvec[2] = 10000;

            Box box(IntVect(0, 0, 0), IntVect(8, 8, 8));
            FArrayBox fab_in(box, 3);
            FArrayBox fab_out(box, 3);
            auto box_pointers = BoxPointers{fab_in, fab_out};
            Cell<double> current_cell(mv_intvec, box_pointers);
            Coordinates<double> coords_metricvars(current_cell, dx,
                                                  kerr_params.center);

            std::cout << "coords: " << coords_metricvars.x << " "
                      << coords_metricvars.y << " " << coords_metricvars.z
                      << std::endl;

            kerr_init.compute_metric_background(metric_vars, coords_metricvars);

            std::cout << "lapse: " << metric_vars.lapse << std::endl;

            std::cout << "shift: " << metric_vars.shift[0] << " "
                      << metric_vars.shift[1] << " " << metric_vars.shift[2]
                      << std::endl;

        } // end of kerr tests
    }

    return 0;
}
