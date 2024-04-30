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
#define ERR 1e-9

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
#include "FourthOrderDerivatives.hpp"


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
#include "MatterEvolutionTest.hpp"
#include "ProcaField.hpp"
#include "KerrdeSitter.hpp"

#include "BaseProcaField.hpp"
#include "ProcaField.hpp"
#include "ProcaFieldTest.hpp"


// c++ includes
#include <algorithm>
#include <iostream>
#include <string>
#include <iomanip>

std::string tab{"\t"}; 
static const auto precision  = std::setprecision(16);
static const int printwidth = 25;




int test1()
{
    std::cout << "\n\n KerrdeSitter test:" << std::endl;

    int failed = 0;

    double dx{0.125};

    // initialize an IntVect object
    IntVect intvect;
    intvect[0] = 1;
    intvect[1] = 0.5;
    intvect[2] = 0.5; 

    // create a coordinate object
    // We set the coordinates to some arbitrary value, but which
    // correspond to the
    //       values in the Mathematica notebook
    Coordinates<double> coords(
        intvect, dx,
        {8,8,8}); // should correspond to the point (x,y,z) = (1,1,1)
    const double x = coords.x;
    const double y = coords.y;
    const double z = coords.z;
    std::cout << tab << "(x,y,z) = (" << x << "," << y << "," << z
                << ")" << std::endl;

    // initialize the kerr params and the kerr object
    KerrdeSitter::params_t kerr_params;
    kerr_params.mass = 1.;
    kerr_params.spin = 0.5;
    kerr_params.cosmo_constant = 0.001;
    kerr_params.center = {8.,8.,8.};
    KerrdeSitter kerr_init(kerr_params, dx);
    std::cout << tab << precision << " kerr mass = " << kerr_params.mass
                << std::endl;
    std::cout << tab  << " kerr spin = " << kerr_params.spin
                << std::endl;
    std::cout << tab
                << " kerr cosmo_constant = " << kerr_params.cosmo_constant
                << std::endl;

    // Discriminant test
    double Q{KerrdeSitter::discriminant<double>(kerr_params)};
    std::cout << tab << "Q = " << Q << std::endl;
    double mathematica_result_Q{0.0175359611526745};
    if (abs(Q - mathematica_result_Q) > ERR)
    {
        std::cout << tab << "Discriminant test failed" << std::endl;
        failed = 1;

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
    mv_intvec[0] = 1;
    mv_intvec[1] = 0.5;
    mv_intvec[2] = 0.5;

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

    auto gamma_UU = TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
    auto chris_phys = TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL;

    std::cout << "lapse: \n" << metric_vars.lapse << std::endl;

    std::cout << "\nshift: \n" << metric_vars.shift[0] << " "
                << metric_vars.shift[1] << " " << metric_vars.shift[2]
                << std::endl;

    std::cout << "\ngamma: \n" << metric_vars.gamma[0][0] << " "
                << metric_vars.gamma[0][1] << " "
                << metric_vars.gamma[0][2] << "\n "
                << metric_vars.gamma[1][0] << " "
                << metric_vars.gamma[1][1] << " "
                << metric_vars.gamma[1][2] << "\n "
                << metric_vars.gamma[2][0] << " "
                << metric_vars.gamma[2][1] << " "
                << metric_vars.gamma[2][2] << std::endl;
    
    std::cout << "\ngammaUU: \n" << gamma_UU[0][0] << " "
                << gamma_UU[0][1] << " "
                << gamma_UU[0][2] << "\n "
                << gamma_UU[1][0] << " "
                << gamma_UU[1][1] << " "
                << gamma_UU[1][2] << "\n "
                << gamma_UU[2][0] << " "
                << gamma_UU[2][1] << " "
                << gamma_UU[2][2] << std::endl;

    std::cout << "\nchris phys: " << std::endl;
    FOR2(i, j)
    {
        FOR1(k)
        {
            std::cout << "  chris_phys[" << i << "][" << j << "][" << k
                        << "] = " << chris_phys[i][j][k];
        }
        std::cout << std::endl;
    }

    std::cout << "\nmetric_vars.d1_lapse: \n" << metric_vars.d1_lapse[0]
                << " " << metric_vars.d1_lapse[1] << " "
                << metric_vars.d1_lapse[2] << std::endl;

    std::cout << "\nmetric_vars.d1_shift: " << std::endl;
    FOR1(i)
    {
        FOR1(j)
        {
            std::cout << "  d1_shift[" << i << "][" << j << "] = " << metric_vars.d1_shift[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nmetric_vars.d1_gamma: \n" << std::endl;
    FOR2(i, j)
    {
        FOR1(k)
        {
            std::cout << "  d1_gamma[" << i << "][" << j << "][" << k
                        << "] = " << metric_vars.d1_gamma[i][j][k];
        }
        std::cout << std::endl;
    }

    std::cout << "\nmetric_vars.k_tensor: " << std::endl;
    FOR2(i, j)
    {
        std::cout << "  k_tensor[" << i << "][" << j
                    << "]= " << metric_vars.K_tensor[i][j];
        std::cout << std::endl;
    }

    std::cout << "\nmetric_vars.K: " << metric_vars.K << std::endl;




    // Verify time derivative of chi vanishes
    std::vector<double> detgamma_deriv { 0,0,0};
    double det_gamma = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
    FOR2(i,j)
    {
        detgamma_deriv[j] += 2 * det_gamma * chris_phys[i][i][j];
    }

    double chi = pow(det_gamma, -1.0/3.0);
    std::vector<double> chi_deriv {0,0,0};
    FOR1(i)
    {
        chi_deriv[i] += detgamma_deriv[i] * (-1.0/3.0) * pow(det_gamma,-4./3.);
    }

    std::cout << "\nchi: " << chi << std::endl;
    std::cout << "\nchi_deriv: " << chi_deriv[0] << " " << chi_deriv[1] << " " << chi_deriv[2] << std::endl;

    double chi_dot = 0;
    FOR1(i)
    {
        chi_dot += metric_vars.shift[i] * chi_deriv[i] - 2./3. * chi * metric_vars.d1_shift[i][i];
    }
    chi_dot += 2./3. * chi * metric_vars.lapse * metric_vars.K;
    std::cout << "\nchi_dot: " << chi_dot << std::endl;

    return 0;

} // end of kerr tests


// This struct should be a copy of the background param_t struct, with specified
// values
struct
{
    double mass = 1.0;
    std::array<double, CH_SPACEDIM> center; // automatically filled in during test2 eval
    double cosmo_constant = 0.001;
    double spin = 0.5;
    double r_plus = 1.868535723767148;

} BackgroundParams;



int test2(bool debug_plots_on = false)
{

    #ifdef _OPENMP
    std::cout << " Number of threads: " << omp_get_max_threads() << std::endl;
    #endif

    int failed{0};                    // flag for failed
    const std::vector<int> resolutions {96, 192, 384};   //resolutions to run at


    const int num_resolutions = 3; //how many of the resolutions to actually use

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
        const int N_GRID{ resolutions[ires] }; // number of cells on each side of box

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

        // create background
        typename KerrdeSitter::params_t bg_params;
        bg_params.mass = BackgroundParams.mass;
        bg_params.spin = BackgroundParams.spin;
        bg_params.cosmo_constant = BackgroundParams.cosmo_constant;
        bg_params.center = center_vector;
        bg_params.r_plus = BackgroundParams.r_plus;
        std::cout << "Mass: " << bg_params.mass << std::endl;
        std::cout << "Spin: " << bg_params.spin << std::endl;
        std::cout << "cosmo_constant: " << bg_params.cosmo_constant << std::endl;
        std::cout << "outer horizon: " << bg_params.r_plus << std::endl;
        std::cout << "Center: " << bg_params.center[0] << " "
                  << bg_params.center[1] << " " << bg_params.center[2]
                  << std::endl;
        std::cout << "dx: " << dx << std::endl;
        KerrdeSitter background_init(bg_params, dx);

        std::cout << "Computing fixed background..." << std::endl;
        // assign background variables to grid
        BoxLoops::loop(AssignFixedBGtoBSSNVars<KerrdeSitter>(background_init,
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
            matter, ccz4_params, dx, sigma, CCZ4RHS<>::USE_BSSN, G_Newt, BackgroundParams.cosmo_constant);

        std::cout << "Numerically computing rhs..." << std::endl;

        // numerically compute RHS
        BoxLoops::loop(matter_ccz4_rhs, fixedbg_FAB, rhs_FAB, disable_simd());

        std::cout << "Analytic computing rhs..." << std::endl;

        // setup analytic Proca field, for computing analytic derivatives
        ProcaField::params_t proca_params = {1, 1, 1};
        ProcaField analytic_matter(background_init, proca_params);

        // compute RHS using analytic expressions
        MatterEvolution<ProcaField, KerrdeSitter> my_an_evolution(
            analytic_matter, background_init, sigma, dx, center_vector);
        BoxLoops::loop(my_an_evolution, fixedbg_FAB, fixedbg_rhs_FAB);

        std::cout << "Excising..." << std::endl;

        // excise the center where values are always large
        ExcisionTest<BaseProcaField<KerrdeSitter, ProcaField>, KerrdeSitter>
            excision(dx, center_vector, background_init);
        BoxLoops::loop(excision, rhs_FAB, rhs_FAB, disable_simd());
        BoxLoops::loop(excision, fixedbg_rhs_FAB, fixedbg_rhs_FAB, disable_simd());

        if (debug_plots_on)
        {
            { //Output parts of the rhs_FAB after subtraction

                std::cout << "In debugging block" << std::endl;
                std::string filename { "output_rhs_res" + std::to_string(ires) + ".txt" };
                std::ofstream outfile;
                outfile.clear();
                outfile.open(filename);
                outfile <<  "#" << std::setw(printwidth) << "x" << std::setw(printwidth) << "y" << std::setw(printwidth) << "z";
                outfile << std::setw(printwidth) <<  "chi";
                outfile << std::setw(printwidth) << "K";
                outfile << std::setw(printwidth) << "lapse" << std::setw(printwidth) << "Avec1" << std::setw(printwidth) << "Avec2" << std::setw(printwidth) << "Avec3" << std::setw(printwidth) << "Ham";
                outfile << std::setw(printwidth) << "Mom1" << std::setw(printwidth) << "Mom2" << std::setw(printwidth) << "Mom3" <<std::endl;
                BoxIterator bit(box);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    IntVect iv = bit();
                    if (iv[1] == N_GRID / 2 && iv[2] == N_GRID / 2)
                    {
                        double x { dx * (iv[0] + 0.5) };
                        double y { dx * (iv[1] + 0.5) };
                        double z { dx * (iv[2] + 0.5) };
                        double out7 = rhs_FAB(iv, c_chi);
                        double out1 = rhs_FAB(iv, c_lapse);
                        double out2 = rhs_FAB(iv, c_Avec1);
                        double out8 = rhs_FAB(iv, c_Avec2);
                        double out9 = rhs_FAB(iv, c_Avec3);
                        double out3 = rhs_FAB(iv, c_Ham);
                        double out4 = rhs_FAB(iv, c_Mom1);
                        double out5 = rhs_FAB(iv, c_Mom2);
                        double out6 = rhs_FAB(iv, c_Mom3);
                        double out10 = rhs_FAB(iv, c_K);
                        
                        outfile << precision << std::setw(printwidth) << x << std::setw(printwidth) << y << std::setw(printwidth) << z ;
                        outfile << std::setw(printwidth) << out7;
                        outfile << std::setw(printwidth) << out10;
                        outfile << std::setw(printwidth) << out1 << std::setw(printwidth) << out2 << std::setw(printwidth)<<out8<< std::setw(printwidth)<<out9 << std::setw(printwidth);
                        outfile << out3 << std::setw(printwidth) << out4 << std::setw(printwidth) << out5;
                        outfile << std::setw(printwidth) << out6 << std::endl;
                    }
                }
                outfile.close();
            }

            { //Output parts of the fixedbg_rhs_FAB after subtraction

                std::cout << "In debugging block" << std::endl;
                std::string filename { "output_fixedbg_rhs_res"+    std::to_string(ires) + ".txt" };
                std::ofstream outfile;
                outfile.clear();
                outfile.open(filename);
                outfile <<  "#" << std::setw(printwidth) << "x" << std::setw(printwidth) << "y" << std::setw(printwidth) << "z";
                outfile << std::setw(printwidth) <<  "chi";
                outfile << std::setw(printwidth) << "K";
                outfile << std::setw(printwidth) << "lapse" << std::setw(printwidth) << "Avec1" << std::setw(printwidth) << "Avec2" << std::setw(printwidth) << "Avec3" << std::setw(printwidth) <<"Ham";
                outfile << std::setw(printwidth) << "Mom1" << std::setw(printwidth) << "Mom2" << std::setw(printwidth) << "Mom3" <<std::endl;
                BoxIterator bit(box);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    IntVect iv = bit();
                    if (iv[1] == N_GRID / 2 && iv[2] == N_GRID / 2)
                    {
                        double x { dx * (iv[0] + 0.5) };
                        double y { dx * (iv[1] + 0.5) };
                        double z { dx * (iv[2] + 0.5) };
                        double out7 = fixedbg_rhs_FAB(iv, c_chi);
                        double out1 = fixedbg_rhs_FAB(iv, c_lapse);
                        double out2 = fixedbg_rhs_FAB(iv, c_Avec1);
                        double out8 = fixedbg_rhs_FAB(iv, c_Avec2);
                        double out9 = fixedbg_rhs_FAB(iv, c_Avec3);
                        double out3 = fixedbg_rhs_FAB(iv, c_Ham);
                        double out4 = fixedbg_rhs_FAB(iv, c_Mom1);
                        double out5 = fixedbg_rhs_FAB(iv, c_Mom2);
                        double out6 = fixedbg_rhs_FAB(iv, c_Mom3);
                        double out10 = fixedbg_rhs_FAB(iv, c_K);
                        
                        outfile << precision << std::setw(printwidth) << x << std::setw(printwidth) << y << std::setw(printwidth) << z ;
                        outfile << std::setw(printwidth) << out7;
                        outfile << std::setw(printwidth) << out10;
                        outfile << std::setw(printwidth) << out1 << std::setw(printwidth) << out2 <<std::setw(printwidth)<<out8<< std::setw(printwidth)<<out9<< std::setw(printwidth);
                        outfile << out3 << std::setw(printwidth) << out4 << std::setw(printwidth) << out5;
                        outfile << std::setw(printwidth) << out6 << std::endl;
                    }
                }
                outfile.close();
            }
        }

        // take the difference between the numerically computed RHS and the
        // analytically computed RHS. should converge to zero for increasing
        // resolution
        rhs_FAB -= fixedbg_rhs_FAB;

        if (debug_plots_on)
        {
            { //Output parts of the rhs_FAB after subtraction

                std::cout << "In debugging block" << std::endl;
                std::string filename { "output_diff_res" + std::to_string(ires) + ".txt" };
                std::ofstream outfile;
                outfile.clear();
                outfile.open(filename);
                outfile <<  "#" << std::setw(printwidth) << "x" << std::setw(printwidth) << "y" << std::setw(printwidth) << "z";
                outfile << std::setw(printwidth) <<  "chi";
                outfile << std::setw(printwidth) << "K";
                outfile << std::setw(printwidth) << "lapse" << std::setw(printwidth) << "Avec1" << std::setw(printwidth) << "Avec2" << std::setw(printwidth) << "Avec3" << std::setw(printwidth) << "Ham";
                outfile << std::setw(printwidth) << "Mom1" << std::setw(printwidth) << "Mom2" << std::setw(printwidth) << "Mom3" <<std::endl;
                BoxIterator bit(box);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    IntVect iv = bit();
                    if (iv[1] == N_GRID / 2 && iv[2] == N_GRID / 2)
                    {
                        double x { dx * (iv[0] + 0.5) };
                        double y { dx * (iv[1] + 0.5) };
                        double z { dx * (iv[2] + 0.5) };
                        double out7 = rhs_FAB(iv, c_chi);
                        double out1 = rhs_FAB(iv, c_lapse);
                        double out2 = rhs_FAB(iv, c_Avec1);
                        double out8 = rhs_FAB(iv, c_Avec2);
                        double out9 = rhs_FAB(iv, c_Avec3);
                        double out3 = rhs_FAB(iv, c_Ham);
                        double out4 = rhs_FAB(iv, c_Mom1);
                        double out5 = rhs_FAB(iv, c_Mom2);
                        double out6 = rhs_FAB(iv, c_Mom3);
                        double out10 = rhs_FAB(iv, c_K);
                        
                        outfile << precision << std::setw(printwidth) << x << std::setw(printwidth) << y << std::setw(printwidth) << z ;
                        outfile << std::setw(printwidth) << out7;
                        outfile << std::setw(printwidth) << out10;
                        outfile << std::setw(printwidth) << out1 << std::setw(printwidth) << out2 << std::setw(printwidth)<<out8<< std::setw(printwidth)<<out9 << std::setw(printwidth);
                        outfile << out3 << std::setw(printwidth) << out4 << std::setw(printwidth) << out5;
                        outfile << std::setw(printwidth) << out6 << std::endl;
                    }
                }
                outfile.close();
            }
        }

        

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

        //check the RHS of the metric variables is zero
        for ( int i { c_chi}; i < c_lapse; i++)
        {
            double max_err = rhs_FAB.norm(max_norm, i, num_comps);
            if (max_err > error_limit)
            {
                std::cout << "RHS FOR COMPONENT "
                          << UserVariables::variable_names[i]
                          << " IS NON ZERO: MAX ERROR = " << max_err
                          << std::endl;
                failed = -1;
            }
            // save the L1 norm for the convergence check - gives sum of abs
            // values
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
    std::vector<std::string> argvs;
    for (int i = 0; i < argc; i++)
    {
        argvs.push_back(argv[i]);
    }

    bool run_test1 = (std::find(argvs.begin(), argvs.end(), "test1") != argvs.end());
    bool run_test2 = (std::find(argvs.begin(), argvs.end(), "test2") != argvs.end());
    bool debug = (std::find(argvs.begin(), argvs.end(), "debug") != argvs.end());



    std::cout << std::string(10, '#') << std::endl;
    std::cout << "Running GRBoondi tests" << std::endl;
    std::cout << std::string(10, '#') << std::endl;

    bool test1res = 0;
    bool test2res = 0;

    if (run_test1)
    {
        std::cout << std::string(10, '#') << std::endl;
        std::cout << "Running GRBoondi tests" << std::endl;
        std::cout << "current tests: KerrdeSitter Simple Test" << std::endl;
        std::cout << std::string(10, '#') << std::endl;
        test1res = test1();
        // Here we test the outputs of the background functions
    }

    if (run_test2)
    {
        std::cout << "\n\n" << std::endl;
        std::cout << std::string(10, '#') << std::endl;
        std::cout << std::string(10, '#') << std::endl;
        std::cout << "current tests: KerrdeSitter Extended Test" << std::endl;
        std::cout << std::string(10, '#') << std::endl;
        std::cout << std::string(10, '#') << std::endl;
        test2res = test2(debug);
    }


    return test2res || test1res;
}
