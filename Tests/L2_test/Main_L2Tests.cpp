

#ifdef _OPENMP
#include <omp.h>
#endif

// background includes
#include "ADMFixedBGVars.hpp"
#include "KerrSchildNew.hpp"
#include "MatterEvolution.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components

// GRChombo includes
#include "BoxLoops.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "SetValue.hpp"

// tests includes
#include "FixedBGProcaField.hpp"
#include "Potential.hpp"

// L2 includes
#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "ProcaField.hpp"

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

// generating random numbers
#include <random>

// c++ includes
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

int test1()
{
#ifdef _OPENMP
    std::cout << " Number of threads: " << omp_get_max_threads() << std::endl;
#endif

    int failed{0};     // flag for failed
    int resolution{0}; // exponent for factor of 2 on default resolution
    int N_GRID{128 *
               (int)pow(2, resolution)}; // number of cells on each side of box

    std::array<double, NUM_VARS> error_norms;

    std::cout << "Creating grid (" << N_GRID << " x " << N_GRID << " x "
              << N_GRID << ")" << std::endl;

    // setup boxes and ghosted boxes
    Box box(IntVect(0, 0, 0), IntVect(N_GRID - 1, N_GRID - 1,
                                      N_GRID - 1)); // The computational box
    Box ghosted_box(
        IntVect(-3, -3, -3),
        IntVect(N_GRID + 2, N_GRID + 2,
                N_GRID + 2)); // The computational box with added ghost cells
    Box double_ghosted_box(
        IntVect(-6, -6, -6),
        IntVect(N_GRID + 5, N_GRID + 5,
                N_GRID +
                    5)); // The computational box with added double ghost cells

    // setup containers for data
    /*     FArrayBox fixedbg_FAB(double_ghosted_box, NUM_VARS); //NUM_Vars from
       UserVariables FArrayBox deriv_fixedbg_FAB(ghosted_box, NUM_VARS);  */
    FArrayBox ref_FAB(double_ghosted_box, NUM_VARS);
    FArrayBox my_FAB(double_ghosted_box, NUM_VARS);
    FArrayBox ref_rhs_FAB(box, NUM_VARS);
    FArrayBox my_rhs_FAB(box, NUM_VARS);

    // set everything to zero
    BoxLoops::loop(SetValue(0.0), ref_FAB, ref_FAB);
    BoxLoops::loop(SetValue(0.0), my_FAB, my_FAB);
    BoxLoops::loop(SetValue(0.0), ref_rhs_FAB, ref_rhs_FAB);
    BoxLoops::loop(SetValue(0.0), my_rhs_FAB, my_rhs_FAB);

    // setup various aspects of the grid
    const double length = 16.0; // coordinate size of the computational box
    const double dx = length / (N_GRID); // grid spacing
    const double center = length / 2.0;  // center of the computational box
    const std::array<double, CH_SPACEDIM> center_vector = {
        center, center, center}; // center of the computational box

    std::cout << "Length: " << length << " Grid spacing: " << dx
              << " Center: " << center << std::endl;

    // create the background
    //  Use a KerrSchildNew background since many of the metric derivatives are
    //  non-zero, testing all terms in the EOM's
    KerrSchildNew::params_t kerr_params{1.0, center_vector, 0.5};
    KerrSchildNew kerr_init(kerr_params, dx);

    // setup the reference class
    Potential::params_t potparams{1.0, 0.0};
    Potential potential(potparams);
    FixedBGProcaField<Potential> test_proca_field(potential, 1.0);
    std::cout << "Test Proca mass: " << potparams.mass << std::endl;
    std::cout << "Test Proca coupling: " << potparams.self_interaction
              << std::endl;

    // setup the GRBoondi L2 generalized Proca class
    ProcaField::params_t myparams{1.0, 1.0, 1.0}; // mass, alpha2,
                                                  // vector_damping
    ProcaField my_proca_field(kerr_init, myparams);

    // Now compute the RHS

    // GRBoondi evolution
    MatterEvolution<ProcaField, KerrSchildNew> my_matter(
        my_proca_field, kerr_init, 0.0, dx,
        center_vector); // set Kreiss-Oliger to zero

    // reference evolution
    MatterEvolution<FixedBGProcaField<Potential>, KerrSchildNew> ref_matter(
        test_proca_field, kerr_init, 0.0, dx,
        center_vector); // set Kreiss-Oliger to zero

    // Now loop over the box and compute the RHS
    std::cout << "Computing GRBoondi's rhs..." << std::endl;
    BoxLoops::loop(my_matter, my_FAB, my_rhs_FAB);
    std::cout << "Computing reference rhs..." << std::endl;
    BoxLoops::loop(ref_matter, ref_FAB, ref_rhs_FAB);

    // now subtract the two
    std::cout << "Computing difference..." << std::endl;
    ref_rhs_FAB -= my_rhs_FAB;

    // checking the results
    const int max_norm = 0;
    const int L1_norm = 1;
    const int num_comps = 1;
    const double error_limit = 1e-10;

    for (int i{c_phi}; i <= c_Z; i++)
    {
        double max_err = ref_rhs_FAB.norm(max_norm, i, num_comps);
        if (max_err > error_limit)
        {
            std::cout << "DERIVATIVES DO NOT MATCH FOR VARIABLE "
                      << UserVariables::variable_names[i] << std::endl;
            failed = -1;
        }
        error_norms[i] = max_err;
    }

    std::cout << " \tMax phi error: " << error_norms[c_phi] << std::endl;
    std::cout << " \tMax Avec1 error: " << error_norms[c_Avec1] << std::endl;
    std::cout << " \tMax Avec2 error: " << error_norms[c_Avec2] << std::endl;
    std::cout << " \tMax Avec3 error: " << error_norms[c_Avec3] << std::endl;
    std::cout << " \tMax Evec1 error: " << error_norms[c_Evec1] << std::endl;
    std::cout << " \tMax Evec2 error: " << error_norms[c_Evec2] << std::endl;
    std::cout << " \tMax Evec3 error: " << error_norms[c_Evec3] << std::endl;
    std::cout << " \tMax Z error: " << error_norms[c_Z] << std::endl;

    return failed;
}

int main(int argc, char *argv[])
{
    std::cout << std::string(10, '#') << std::endl;
    std::cout << "Running GRBoondi tests" << std::endl;
    std::cout << "current tests: L2 tests" << std::endl;
    std::cout << std::string(10, '#') << std::endl;

    // Here we test the L2 routines
    int res = test1();

    if (res == 0)
    {
        std::cout << "GRBoondi's RHS matches reference ..... PASS" << std::endl;
    }
    else
    {
        std::cout << "GRBoondi's RHS does not match reference ..... FAILED"
                  << std::endl;
    }

    return res;
}
