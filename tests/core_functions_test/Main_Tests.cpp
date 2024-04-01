
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
            std::vector<std::vector<double>> adxx_sol = {
                {0.248963113140091, 0.136392981253951, -0.110832940095943},
                {0.136392981253951, 0.248963113140091, -0.110832940095943},
                {-0.110832940095943, -0.110832940095943,  -0.265837195347114}
            };

            std::vector<std::vector<std::vector<double>>> betadxx_sol = {
                {
                    {-0.918081421262146,-0.437063708465753,0.358635632821701},
                    {-0.437063708465753,-1.31597450368002,0.118126776423505},
                    {0.358635632821702, 0.118126776423505, 0.699284505275556}
                },
                {
                    {0.701676194416934, -0.177234600797329, 0.131937016697419},
                    {-0.177234600797329, -0.658252313593723, -1.01803328600267},
                    {0.131937016697419, -1.01803328600267, -0.592906114534130}
                },
                {
                    {-0.217600047027972, 0.488948715795823, 0.435918368878509},
                    {0.488948715795823, -0.217600047027972, 0.435918368878509},
                    {0.435918368878509, 0.435918368878509,  -1.20065844107154}
                }
            };

            double err { 1e-14 };

            bool alphapass { true };
            bool shiftpass { true };

            //output
            std::cout << std::setprecision (15) << std::endl;
            FOR2(i,j)
            {

                if (std::abs(metric_vars.d2_lapse[i][j] - adxx_sol[i][j]) > err)
                {
                    std::cout << tab << "d2_lapse_d" << i << "_d" << j << " = " << metric_vars.d2_lapse[i][j] << " != " << adxx_sol[i][j] << std::endl;
                    alphapass = false;
                }

                FOR1(k)
                {
                    if (std::abs(metric_vars.d2_shift[i][j][k] - betadxx_sol[i][j][k]) > err)
                    {
                        std::cout << tab << "d2_shift^" << i << "_d" << j << "_d" << k << " = " << metric_vars.d2_shift[i][j][k] << " != " << betadxx_sol[i][j][k] << std::endl;
                        shiftpass = false;
                    }
                }
            }

            if (alphapass)
            {
                std::cout << tab << "d2_lapse test passed" << std::endl;
            } else {
                std::cout << tab << "d2_lapse test failed" << std::endl;
            }

            if (shiftpass)
            {
                std::cout << tab << "d2_shift test passed" << std::endl;
            } else {
                std::cout << tab << "d2_shift test failed" << std::endl;
            }

        }//end of kerr tests

    }



    return 0;
}