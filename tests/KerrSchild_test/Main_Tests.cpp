
// background includes
#include "KerrSchild.hpp"
#include "ADMFixedBGVars.hpp"
#include "KerrSchildFixedBG.hpp"

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
    std::cout << "Running GRBoondi tests" << std::endl;
    std::cout << "current tests: KerrSchild" << std::endl;
    std::cout << std::string(10, '#') << std::endl;

    std::string tab { "\t" };

    //Here we test the outputs of the background functions
    {

        //KerrSchild
        {
            std::cout << "\n\n KerrSchild test:" << std::endl;

            double dx { 1e-3 };

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
            //init reference class
            KerrSchildFixedBG::params_t ref_params = {kerr_params.mass, kerr_params.center, kerr_params.spin};
            KerrSchildFixedBG kerr_ref (ref_params, dx);

            //Excision test
            {
               //setup container objects for testing excision code
                std::vector<double> ref_excision;
                std::vector<double> my_excision;
                IntVect excision_intvect(1.,1.,1.);

                //First check BH params match
                if (kerr_init.m_params.mass != kerr_ref.m_params.mass)
                {
                    std::cout << "Mass test failed" << std::endl;
                }
                if (kerr_init.m_params.spin != kerr_ref.m_params.spin)
                {
                    std::cout << "Spin test failed" << std::endl;
                }

                bool ExcisionPass { true };
                double buffer {0.9554};
                for (int i {0}; i<3/dx; i++)
                {
                    excision_intvect[0] = 1.0;
                    excision_intvect[1] = 1.0;
                    excision_intvect[2] += i;

                    Box box(IntVect(0,0,0), IntVect(8,8,8));
                    FArrayBox fab_in(box,3);
                    FArrayBox fab_out(box,3);
                    auto box_pointers = BoxPointers{fab_in, fab_out};
                    Cell<double> current_cell(excision_intvect, box_pointers);
                    const Coordinates<double> coords(current_cell, dx, kerr_params.center);
                    bool kerr_ref_excise_result { kerr_ref.excise(current_cell)<buffer };
                    bool kerr_init_excise_result { kerr_init.check_if_excised(coords,buffer)  };
                    if (kerr_ref_excise_result != kerr_init_excise_result)
                    {
                        
                        std::cout << "Excision test failed" << std::endl;
                        std::cout << "Coordinate location: " << coords << std::endl;
                        std::cout << "kerr reference output: " << kerr_ref_excise_result << std::endl;
                        std::cout << "kerr my output: " << kerr_init_excise_result << std::endl;
                        ExcisionPass = false;
                    }

                }
                std::cout << "Excision result: " << ExcisionPass << std::endl;
            }
            

           

        }//end of kerr tests

    }



    return 0;
}
