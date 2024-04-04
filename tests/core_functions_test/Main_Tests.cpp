
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
            Coordinates<double> coords (intvect, dx, {0,0,0}); // should correspond to the point (x,y,z) = (1,1,1)
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

            //These solutions come directly from mathematica
	    double adxx_sol_sum { 0.06154323305720 };
	    double betadxx_sol_sum { -3.0458056200456 };
	    double gammadxx_sol_sum { -22.03836675054085 };
            double err { 1e-10 };

            bool alphapass { true };
            bool shiftpass { true };
	    bool gammapass { true };

            //output
            std::cout << std::setprecision (15) << std::endl;
	    double adxx { 0. };
	    double betadxx { 0. };
	    double gammadxx { 0. };
            FOR2(i,j)
            {
		adxx += metric_vars.d2_lapse[i][j];
                FOR1(k)
                {
		    betadxx += metric_vars.d2_shift[i][j][k];
		    FOR1(l)
		    {
			gammadxx += metric_vars.d2_gamma[i][j][k][l];
		    }
                }
            }
	    if (adxx - adxx_sol_sum > err)
	    {
		alphapass = false;
	    } 
	    if (betadxx - betadxx_sol_sum > err)
	    {
		shiftpass = false;
	    };
	    if (gammadxx - gammadxx_sol_sum > err)
	    {
		gammapass = false;
	    }

            if (alphapass)
            {
                std::cout << tab << "d2_lapse test passed" << std::endl;
            } else {
                std::cout << tab << "d2_lapse test failed" << std::endl;
		std::cout << tab << "adxx_sum = " << adxx << std::endl;
            }

            if (shiftpass)
            {
                std::cout << tab << "d2_shift test passed" << std::endl;
            } else {
                std::cout << tab << "d2_shift test failed" << std::endl;
		std::cout << tab << "betadxx_sum = " << betadxx << std::endl;
            }

	    if (gammapass)
	    {
		std::cout << tab << "d2_gamma test passed" << std::endl;
	    } else {
		std::cout << tab << "d2_gamma test failed" << std::endl;
		std::cout << tab << "gammadxx_sum = " << gammadxx << std::endl;
	    };

        }//end of kerr tests

    }



    return 0;
}
