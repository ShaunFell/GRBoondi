/*
GRBoondi

Some of these tests were borrowed from the GRDzhadzha test suite. Please see
https://github.com/GRTLCollaboration/GRDzhadzha/tree/main/Tests/KerrBHScalarTest
for more details
*/

#ifdef _OPENMP
#include <omp.h>
#endif

#define NUM_ERR_LIM 0.1

// misc includes
#include "AssignFixedBGtoBSSNVars.hpp"
#include "ExcisionTest.hpp"

// GRChombo includes
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "GammaCalculator.hpp"
#include "MatterCCZ4.hpp"
#include "MatterCCZ4RHS.hpp"
#include "NewConstraints.hpp"

// Chombo includes
#include "BoxLoops.hpp"
#include "Cell.hpp"
#include "ComputePack.hpp"
#include "DebuggingTools.hpp"
#include "FArrayBox.H"
#include "IntVect.H"
#include "SetValue.hpp"
#include "UsingNamespace.H"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/time.h>

// background includes
#include "ADMFixedBGVars.hpp"
#include "BaseProcaFieldTest.hpp"
#include "KerrSchild.hpp"
#include "KerrSchildFixedBG.hpp"
#include "MatterEvolution.hpp"
#include "ProcaField.hpp"

#include "BaseProcaField.hpp"
#include "ProcaField.hpp"
#include "ProcaFieldTest.hpp"

// c++ includes
#include <algorithm>
#include <iostream>
#include <string>

std::string tab{"\t"};

void test1()
{
    std::cout << "\n\n KerrSchild test:" << std::endl;

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
    std::cout << tab << "(x,y,z) = (" << x << "," << y << "," << z << ")"
              << std::endl;

    // initialize the kerr params and the kerr object
    KerrSchild::params_t kerr_params;
    kerr_params.mass = 1.0;
    kerr_params.spin = 0.5;
    // kerr_params.center = {0.0, 0.0, 0.0};
    KerrSchild kerr_init(kerr_params, dx);
    std::cout << tab << " kerr mass = " << kerr_params.mass << std::endl;
    std::cout << tab << " kerr spin = " << kerr_params.spin << std::endl;

    // instantiate the metric variables
    ADMFixedBGVars::Vars<double> metric_vars;

    // the real work
    // init reference class
    KerrSchildFixedBG::params_t ref_params = {
        kerr_params.mass, kerr_params.center, kerr_params.spin};
    KerrSchildFixedBG kerr_ref(ref_params, dx);

    // Excision test
    {
        // setup container objects for testing excision code
        std::vector<double> ref_excision;
        std::vector<double> my_excision;
        IntVect excision_intvect(1., 1., 1.);

        // First check BH params match
        if (kerr_init.m_params.mass != kerr_ref.m_params.mass)
        {
            std::cout << "Mass test failed" << std::endl;
        }
        if (kerr_init.m_params.spin != kerr_ref.m_params.spin)
        {
            std::cout << "Spin test failed" << std::endl;
        }

        bool ExcisionPass{true};
        double buffer{0.9554};
        for (int i{0}; i < 3 / dx; i++)
        {
            excision_intvect[0] = 1.0;
            excision_intvect[1] = 1.0;
            excision_intvect[2] += i;

            Box box(IntVect(0, 0, 0), IntVect(8, 8, 8));
            FArrayBox fab_in(box, 3);
            FArrayBox fab_out(box, 3);
            auto box_pointers = BoxPointers{fab_in, fab_out};
            Cell<double> current_cell(excision_intvect, box_pointers);
            const Coordinates<double> coords(current_cell, dx,
                                             kerr_params.center);
            bool kerr_ref_excise_result{kerr_ref.excise(current_cell) < buffer};
            bool kerr_init_excise_result{
                kerr_init.check_if_excised(coords, buffer)};
            if (kerr_ref_excise_result != kerr_init_excise_result)
            {

                std::cout << "Excision test failed" << std::endl;
                std::cout << "Coordinate location: " << coords << std::endl;
                std::cout << "kerr reference output: " << kerr_ref_excise_result
                          << std::endl;
                std::cout << "kerr my output: " << kerr_init_excise_result
                          << std::endl;
                ExcisionPass = false;
            }
        }
        std::cout << "Excision passed: " << std::boolalpha << ExcisionPass
                  << std::endl
                  << std::endl
                  << std::endl;
    }
};

// This struct should be a copy of the background param_t struct, with specified
// values
struct
{
    double mass = 1.0;
    std::array<double, CH_SPACEDIM>
        center; // automatically filled in during test2 eval
    double spin = 0.5;
} BackgroundParams;

template <class background_t> int test2()
{

    #ifdef _OPENMP
    std::cout << " Number of threads: " << omp_get_max_threads() << std::endl;
    #endif

    int failed{0};                    // flag for failed
    const bool debug_plots_on{false}; // export data for plotting
    const int num_resolutions{2};     // number of resolutions to run at

    // vector of norms for convergence checking
    std::array<std::array<double, NUM_VARS>, num_resolutions> error_norms;

    for (int ires{0}; ires < num_resolutions; ++ires)
    {
        std::cout << std::endl << std::endl;
        std::cout << "*******************" << std::endl;
        std::cout << "Resolution: " << ires << std::endl;
        std::cout << "*******************" << std::endl;

        // fill current error norm array with zeroes
        error_norms[ires].fill(0.0);

        // setup the array boxes for various inputs and outputs
        const int N_GRID{
            96 * (int)pow(2, ires)}; // number of cells on each side of box

        std::cout << "Creating grid (" << N_GRID << " x " << N_GRID << " x "
                  << N_GRID << ")" << std::endl;
        // setup boxes
        Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1,
                                          N_GRID - 1)); // The computational box
        Box ghosted_box(
            IntVect(-3, -3, -3),
            IntVect(N_GRID + 2, N_GRID + 2,
                    N_GRID +
                        2)); // The computational box with added ghost cells
        Box double_ghosted_box(
            IntVect(-6, -6, -6),
            IntVect(
                N_GRID + 5, N_GRID + 5,
                N_GRID +
                    5)); // The computational box with added double ghost cells

        // setup containers for data
        FArrayBox fixedbg_FAB(double_ghosted_box, NUM_VARS);
        FArrayBox deriv_fixedbg_FAB(ghosted_box, NUM_VARS);
        FArrayBox rhs_FAB(box, NUM_VARS);
        FArrayBox fixedbg_rhs_FAB(box, NUM_VARS);

        // loop over box cells and set everything to zero
        BoxLoops::loop(make_compute_pack(SetValue(0.0)), fixedbg_FAB,
                       fixedbg_FAB);
        BoxLoops::loop(make_compute_pack(SetValue(0.0)), deriv_fixedbg_FAB,
                       deriv_fixedbg_FAB);
        BoxLoops::loop(make_compute_pack(SetValue(0.0)), rhs_FAB, rhs_FAB);
        BoxLoops::loop(make_compute_pack(SetValue(0.0)), fixedbg_rhs_FAB,
                       fixedbg_rhs_FAB);

        // setup various aspects of the grid
        const double length = 16.0; // coordinate size of the computational box
        const double dx = length / (N_GRID); // grid spacing
        const double center = length / 2.0;  // center of the computational box
        const std::array<double, CH_SPACEDIM> center_vector = {
            center, center, center}; // center of the computational box
        BackgroundParams.center = center_vector;

        // create background
        typename background_t::params_t bg_params;
        bg_params.mass = BackgroundParams.mass;
        bg_params.spin = BackgroundParams.spin;
        bg_params.center = center_vector;
        std::cout << "Mass: " << bg_params.mass << std::endl;
        std::cout << "spin: " << bg_params.spin << std::endl;
        std::cout << "Center: " << bg_params.center[0] << " "
                  << bg_params.center[1] << " " << bg_params.center[2]
                  << std::endl;
        background_t background_init(bg_params, dx);

        std::cout << "Computing fixed background..." << std::endl;
        // assign background variables to grid
        BoxLoops::loop(AssignFixedBGtoBSSNVars<background_t>(background_init,
                                                             dx, center_vector),
                       fixedbg_FAB, fixedbg_FAB);
        GammaCalculator gamamcalc(dx);
        BoxLoops::loop(gamamcalc, fixedbg_FAB, deriv_fixedbg_FAB);

        // add derivatives to variables
        fixedbg_FAB += deriv_fixedbg_FAB;

        std::cout << "computing constraints..." << std::endl;

        // compute the hamiltonian and momentum constraints and put them in the
        // RHS fab
        BoxLoops::loop(Constraints(dx, c_Ham, Interval(c_Mom1, c_Mom3)),
                       fixedbg_FAB, rhs_FAB, box);

        // set CCZ4 parameters
        const double G_Newt{0.}; // dont evolve the background
        const double sigma{0.};  // turn off Kreiss-Oliger dissipation
        CCZ4RHS<>::params_t ccz4_params;
        ccz4_params.kappa1 = 0.0;
        ccz4_params.kappa2 = 0.0;
        ccz4_params.kappa3 = 0.0;
        ccz4_params.lapse_coeff = 0.0; // no evolution lapse or shift
        ccz4_params.shift_Gamma_coeff = 0.0;

        // setup matter class
        ProcaFieldTest::params_t proca_params_test = {1, 1, 1};
        ProcaFieldTest matter(background_init, proca_params_test);

        // setup matterccz4 rhs with matter class
        MatterCCZ4RHS<ProcaFieldTest> matter_ccz4_rhs(
            matter, ccz4_params, dx, sigma, CCZ4RHS<>::USE_BSSN, G_Newt);

        std::cout << "Numerically computing rhs..." << std::endl;

        // numerically compute RHS
        BoxLoops::loop(matter_ccz4_rhs, fixedbg_FAB, rhs_FAB);

        std::cout << "Analytic computing rhs..." << std::endl;

        // setup analytic Proca field, for computing analytic derivatives
        ProcaField::params_t proca_params = {1, 1, 1};
        ProcaField analytic_matter(background_init, proca_params);

        // compute RHS using analytic expressions
        MatterEvolution<ProcaField, background_t> my_an_evolution(
            analytic_matter, background_init, sigma, dx, center_vector);
        BoxLoops::loop(my_an_evolution, fixedbg_FAB, fixedbg_rhs_FAB);

        // take the difference between the numerically computed RHS and the
        // analytically computed RHS. should converge to zero for increasing
        // resolution
        rhs_FAB -= fixedbg_rhs_FAB;

        std::cout << "Excising..." << std::endl;

        // excise the center where values are always large
        ExcisionTest<BaseProcaField<background_t, ProcaField>, background_t>
            excision(dx, center_vector, background_init);
        BoxLoops::loop(excision, rhs_FAB, rhs_FAB, disable_simd());

        std::cout << "Checking results..." << std::endl;

        // Now we check the results
        const int max_norm = 0;
        const int L1_norm = 1;
        const int num_comps = 1;
        const double error_limit = NUM_ERR_LIM;

        // check that the initial data satisfies the constraints
        for (int i{c_Ham}; i <= c_Mom3; i++)
        {
            double max_err = rhs_FAB.norm(max_norm, i, num_comps);
            if (max_err > error_limit)
            {
                std::cout << "CONSTRAINT " << UserVariables::variable_names[i]
                          << " IS NON ZERO: MAX ERROR = " << max_err
                          << std::endl;
                failed = -1;
            }

            // save the L1 norm for convergence checking
            error_norms[ires][i] =
                rhs_FAB.norm(L1_norm, i, num_comps) * pow(N_GRID, -3);
        }

        // compare the numerically calculated and the analytically calculated
        // RHS. This checks the expressions for the derivatives of the analytic
        // background variables (d1_gamma, d1_lapse, etc.)
        for (int i{c_phi}; i <= c_Z; i++)
        {
            double max_err = rhs_FAB.norm(max_norm, i, num_comps);
            if (max_err > error_limit)
            {
                std::cout
                    << "ANALYTIC MATTER VARS RHS FOR "
                    << UserVariables::variable_names[i]
                    << " DOES NOT MATCH FINITE DIFFERENCE RHS: MAX ERROR = "
                    << max_err << std::endl;
                failed = -1;
            }

            // save for convergence checking
            error_norms[ires][i] =
                rhs_FAB.norm(L1_norm, i, num_comps) * pow(N_GRID, -3);
        }

    } // end of resolution loop

    // check convergence for increasing resolution
    double min_convergence_factor = 16.0;
    for (int i{0}; i < NUM_VARS; i++)
    {
        for (int ires{0}; ires < num_resolutions - 1; ires++)
        {
            double hi_res_norm = error_norms[ires + 1][i];
            double lo_res_norm = error_norms[ires][i];
            // ignore exact zero values
            if (abs(hi_res_norm) < 1e-16 && abs(lo_res_norm) < 1e-16)
            {
                lo_res_norm = 1e-8;
                hi_res_norm = 1e-10;
            }
            double convergence_factor = lo_res_norm / hi_res_norm;
            if (convergence_factor < min_convergence_factor)
            {
                min_convergence_factor = convergence_factor;
            }
            if (convergence_factor < 11)
            {
                failed = -1;
                std::cout << "CONVERGENCE FACTOR FOR COMPONENT "
                          << UserVariables::variable_names[i] << " ON LEVEL "
                          << ires << " IS LOW: VALUE = " << convergence_factor
                          << " " << hi_res_norm << " " << lo_res_norm
                          << std::endl;
            }
        }
    }

    if (failed == 0)
    {
        std::cout << "The minimum convergence factor was "
                  << min_convergence_factor << std::endl;
        std::cout << "Fixed Background test ........ passed" << std::endl;
    }
    else
    {
        std::cout << "The minimum convergence factor was "
                  << min_convergence_factor << std::endl;
        std::cout << "Fixed Background test ........ failed!!!" << std::endl;
    }

    return failed;
}

int main(int argc, char *argv[])
{
    std::cout << std::string(10, '#') << std::endl;
    std::cout << "Running GRBoondi tests" << std::endl;
    std::cout << std::string(10, '#') << std::endl;

    // Here we test the outputs of the background functions

    std::cout << std::string(10, '#') << std::endl;
    std::cout << std::string(10, '#') << std::endl;
    std::cout << "current tests: KerrSchild" << std::endl;
    std::cout << std::string(10, '#') << std::endl;
    std::cout << std::string(10, '#') << std::endl;
    test1();

    std::cout << std::string(10, '#') << std::endl;
    std::cout << std::string(10, '#') << std::endl;
    std::cout << "current tests: KerrSchild Extended Test" << std::endl;
    std::cout << std::string(10, '#') << std::endl;
    std::cout << std::string(10, '#') << std::endl;
    int test2res = test2<KerrSchild>();

    return test2res;
}
