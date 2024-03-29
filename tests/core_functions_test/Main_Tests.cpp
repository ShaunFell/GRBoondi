
// background includes
#include "KerrSchild.hpp"
#include "ADMFixedBGVars.hpp"

//GRChombo includes
#include "Cell.hpp"
#include "Coordinates.hpp"

//Chombo includes
#include "IntVect.H"
#include "UsingNamespace.H"


//c++ includes
#include <iostream>
#include <string>
#include <algorithm>



int main(int argc, char *argv[])
{
    std::cout << std::string(10, '#') << std::endl;
    std::cout << "Running GRMilijun tests" << std::endl;
    std::cout << "current tests: KerrSchild" << std::endl;
    std::cout << std::string(10, '#') << std::endl;

    std::string tab { "\t" };

    //Here we test the outputs of the background functions
    {

        //KerrSchild
        {
            std::cout << "\n\n KerrSchild test:" << std::endl;

            double dx { 1.0 };

            //initialize an IntVect object
            IntVect intvect;
            intvect[0] = 0.5;
            intvect[1] = 0.5;
            intvect[2] = 0.5;

            //create a coordinate object
            //We set the coordinates to some arbitrary value, but which correspond to the
            //      values in the Mathematica notebook
            Coordinates<double> coords (intvect, dx); // should correspond to the point (x,y,z) = (1,1,1)
            const double x = coords.x;
            const double y = coords.y;
            const double z = coords.z;
            std::cout << tab << "(x,y,z) = (" << x << "," << y << "," << z << ")" << std::endl;


            //initialize the kerr params and the kerr object
            KerrSchild::params_t kerr_params;
            kerr_params.mass = 1.0;
            kerr_params.spin = 0.5;
            //kerr_params.center = {0.0, 0.0, 0.0};
            KerrSchild kerr_init(kerr_params, dx);
            std::cout << tab << " kerr mass = " << kerr_params.mass << std::endl;
            std::cout << tab << " kerr spin = " << kerr_params.spin << std::endl;

            //instantiate the metric variables
            ADMFixedBGVars::Vars<double> metric_vars;

            
            //the real work

            //compute the second derivatives
            kerr_init.compute_2nd_derivatives(metric_vars, coords);

            //output
            FOR2(i,j)
            {
                std::cout << tab << "d2_lapse_d" << i << "_d" << j << " = " << metric_vars.d2_lapse[i][j] << std::endl;

                FOR1(k)
                {
                    std::cout << tab << "d2_shift^" << i << "_d" << j << "_d" << k << " = " << metric_vars.d2_shift[i][j][k] << std::endl;
                }
            }

        }//end of kerr tests

    }



    return 0;
}