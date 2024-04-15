
// background includes
#include "ADMFixedBGVars.hpp"
#include "KerrSchild.hpp"

// GRChombo includes
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "TensorAlgebra.hpp"

// tests includes
#include "FixedBGProcaField.hpp"
#include "Potential.hpp"

// L2 includes
#include "ADMProcaVars.hpp"
#include "BaseProcaField.hpp"
#include "ProcaField.hpp"

// Chombo includes
#include "IntVect.H"
#include "UsingNamespace.H"

// generating random numbers
#include <random>

// c++ includes
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
    std::cout << std::string(10, '#') << std::endl;
    std::cout << "Running GRBoondi tests" << std::endl;
    std::cout << "current tests: L2 tests" << std::endl;
    std::cout << std::string(10, '#') << std::endl;

    std::string tab{"\t"};

    // Here we test the L2 routines
    {
        double dx{1.0};

        std::srand(std::time(nullptr));

        ADMFixedBGVars::Vars<double> metric_vars;
        FixedBGProcaField<Potential>::Vars<double> test_matter_vars;
        FixedBGProcaField<Potential>::Vars<Tensor<1, double>>
            test_matter_diff1_vars;
        FixedBGProcaField<Potential>::Diff2Vars<Tensor<2, double>>
            test_matter_diff2_vars;

        ADMProcaVars::MatterVars<double> my_matter_vars;
        ADMProcaVars::MatterVars<Tensor<1, double>> my_matter_diff1_vars;
        ADMProcaVars::Diff2MatterVars<Tensor<2, double>> my_matter_diff2_vars;

        test_matter_vars.phi = rand() % 10;
        test_matter_vars.Z = rand() % 10;
        my_matter_vars.phi = test_matter_vars.phi;
        my_matter_vars.Z = test_matter_vars.Z;

        metric_vars.lapse = rand() % 10;
        metric_vars.K = rand() % 10;
        FOR1(i)
        {
            metric_vars.shift[i] = rand() % 10;
            metric_vars.d1_lapse[i] = rand() % 10;
            test_matter_vars.Avec[i] = rand() % 10;
            test_matter_vars.Evec[i] = rand() % 10;
            test_matter_diff1_vars.phi[i] = rand() % 10;
            test_matter_diff1_vars.Z[i] = rand() % 10;
            my_matter_vars.Avec[i] = test_matter_vars.Avec[i];
            my_matter_vars.Evec[i] = test_matter_vars.Evec[i];
            my_matter_diff1_vars.phi[i] = test_matter_diff1_vars.phi[i];
            my_matter_diff1_vars.Z[i] = test_matter_diff1_vars.Z[i];

            FOR1(j)
            {
                metric_vars.d1_shift[i][j] = rand() % 10;
                test_matter_diff1_vars.Avec[i][j] = rand() % 10;
                test_matter_diff1_vars.Evec[i][j] = rand() % 10;
                my_matter_diff1_vars.Avec[i][j] =
                    test_matter_diff1_vars.Avec[i][j];
                my_matter_diff1_vars.Evec[i][j] =
                    test_matter_diff1_vars.Evec[i][j];
                if (j >= i)
                {
                    metric_vars.K_tensor[i][j] = rand() % 10;
                    metric_vars.gamma[i][j] = rand() % 10;
                    metric_vars.d2_lapse[i][j] = rand() % 10;
                }
                else
                {
                    metric_vars.K_tensor[i][j] = metric_vars.K_tensor[j][i];
                    metric_vars.gamma[i][j] = metric_vars.gamma[j][i];
                    metric_vars.d2_lapse[i][j] = metric_vars.d2_lapse[j][i];
                }

                FOR1(k)
                {
                    metric_vars.d1_gamma[i][j][k] = rand() % 10;
                    if (k >= j)
                    {
                        metric_vars.d2_shift[i][j][k] = rand() % 10;
                        test_matter_diff2_vars.Avec[i][j][k] = rand() % 10;
                    }
                    else
                    {
                        metric_vars.d2_shift[i][j][k] =
                            metric_vars.d2_shift[i][k][j];
                        test_matter_diff2_vars.Avec[i][j][k] =
                            test_matter_diff2_vars.Avec[i][k][j];
                    }

                    my_matter_diff2_vars.Avec[i][j][k] =
                        test_matter_diff2_vars.Avec[i][j][k];

                    FOR1(m)
                    {
                        if (m >= k)
                        {
                            metric_vars.d2_gamma[i][j][k][m] = rand() % 10;
                        }
                        else
                        {
                            metric_vars.d2_gamma[i][j][k][m] =
                                metric_vars.d2_gamma[i][j][m][k];
                        }
                    }
                }
            }
        }

        const auto gamma_UU =
            TensorAlgebra::compute_inverse_sym(metric_vars.gamma);
        const Tensor<3, double> chris_phys =
            TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU)
                .ULL;

        double proca_mass = 1.0;
        double vec_damping = 1.0;

        Potential::params_t potparams{proca_mass, 0.0};
        Potential potential(potparams);
        FixedBGProcaField<Potential> test_proca_field(potential, vec_damping);

        KerrSchild::params_t kerr_params;
        kerr_params.mass = 1;
        kerr_params.spin = 0.5;
        kerr_params.center = {0.5, 0.5, 0.5};

        KerrSchild kerr_bh{kerr_params, dx};

        emtensor_t<double> Test_EMTest = test_proca_field.compute_emtensor(
            test_matter_vars, metric_vars, test_matter_diff1_vars, gamma_UU,
            chris_phys);

        // My proca field
        ProcaField::params_t myparams{proca_mass, 1.0, vec_damping};
        ProcaField my_proca_field(kerr_bh, myparams);
        emtensor_t<double> my_EMTest = my_proca_field.compute_emtensor(
            my_matter_vars, metric_vars, my_matter_diff1_vars,
            my_matter_diff2_vars, my_matter_vars);
        double testAvecsum{0.};
        double myAvecsum{0.};
        double testEvecsum{0.};
        double myEvecsum{0.};
        double testd1Avecsum{0.};
        double myd1Avecsum{0.};
        double testd1Evecsum{0.};
        double myd1Evecsum{0.};
        double testd2Avecsum{0.};
        double myd2Avecsum{0.};

        FOR1(i)
        {
            testAvecsum += test_matter_vars.Avec[i];
            myAvecsum += my_matter_vars.Avec[i];
            testEvecsum += test_matter_vars.Evec[i];
            myEvecsum += my_matter_vars.Evec[i];
            FOR1(j)
            {
                testd1Avecsum += test_matter_diff1_vars.Avec[i][j];
                myd1Avecsum += my_matter_diff1_vars.Avec[i][j];
                testd1Evecsum += test_matter_diff1_vars.Evec[i][j];
                myd1Evecsum += my_matter_diff1_vars.Evec[i][j];
                FOR1(k)
                {
                    testd2Avecsum += test_matter_diff2_vars.Avec[i][j][k];
                    myd2Avecsum += my_matter_diff2_vars.Avec[i][j][k];
                }
            }
        }

        printf("test phi: %0.8f\n", test_matter_vars.phi);
        printf("my phi: %0.8f\n", my_matter_vars.phi);
        printf("test Z: %0.8f\n", test_matter_vars.Z);
        printf("my Z: %0.8f\n", my_matter_vars.Z);
        printf("test sum Avec: %0.8f\n", testAvecsum);
        printf("my sum Avec: %0.8f\n", myAvecsum);
        printf("test sum Evec: %0.8f\n", testEvecsum);
        printf("my sum Evec: %0.8f\n", myEvecsum);
        printf("test sum d1Evec: %0.8f\n", testd1Evecsum);
        printf("my sum d1Evec: %0.8f\n", myd1Evecsum);
        printf("test sum d1Avec: %0.8f\n", testd1Avecsum);
        printf("my sum d1Avec: %0.8f\n", myd1Avecsum);
        printf("test sum d2Avec: %0.8f\n", testd2Avecsum);
        printf("my sum d2Avec: %0.8f\n", myd2Avecsum);

        std::cout << std::setprecision(16) << std::endl;
        std::cout << "Test rho: " << Test_EMTest.rho << std::endl;
        std::cout << "My rho:   " << my_EMTest.rho << std::endl;
        std::cout << "Test sum Pi: "
                  << Test_EMTest.Si[0] + Test_EMTest.Si[1] + Test_EMTest.Si[2]
                  << std::endl;
        std::cout << "My sum Pi:   "
                  << my_EMTest.Si[0] + my_EMTest.Si[1] + my_EMTest.Si[2]
                  << std::endl;
        std::cout << "Test sum Sij: "
                  << Test_EMTest.Sij[0][0] + Test_EMTest.Sij[0][1] +
                         Test_EMTest.Sij[0][2] + Test_EMTest.Sij[1][0] +
                         Test_EMTest.Sij[1][1] + Test_EMTest.Sij[1][2] +
                         Test_EMTest.Sij[2][1] + Test_EMTest.Sij[2][2] +
                         Test_EMTest.Sij[2][0]
                  << std::endl;
        std::cout << "My sum Sij:   "
                  << my_EMTest.Sij[0][0] + my_EMTest.Sij[0][1] +
                         my_EMTest.Sij[0][2] + my_EMTest.Sij[1][0] +
                         my_EMTest.Sij[1][1] + my_EMTest.Sij[1][2] +
                         my_EMTest.Sij[2][1] + my_EMTest.Sij[2][2] +
                         my_EMTest.Sij[2][0]
                  << std::endl;

        FixedBGProcaField<Potential>::Vars<double> test_RHS;
        ADMProcaVars::MatterVars<double> my_RHS;

        my_proca_field.matter_rhs(my_RHS, my_matter_vars, metric_vars,
                                  my_matter_diff1_vars, my_matter_diff2_vars,
                                  my_matter_vars);
        test_proca_field.matter_rhs(test_RHS, test_matter_vars, metric_vars,
                                    test_matter_diff1_vars,
                                    test_matter_diff2_vars, test_matter_vars);

        std::cout << "Test phi dot: " << test_RHS.phi << std::endl;
        std::cout << "My phi dot:   " << my_RHS.phi << std::endl;
        std::cout << "Test Avec dot sum: "
                  << test_RHS.Avec[0] + test_RHS.Avec[1] + test_RHS.Avec[2]
                  << std::endl;
        std::cout << "My Avec dot sum:   "
                  << my_RHS.Avec[0] + my_RHS.Avec[1] + my_RHS.Avec[2]
                  << std::endl;
        std::cout << "Test Evec dot sum: "
                  << test_RHS.Evec[0] + test_RHS.Evec[1] + test_RHS.Evec[2]
                  << std::endl;
        std::cout << "My Evec dot sum:   "
                  << my_RHS.Evec[0] + my_RHS.Evec[1] + my_RHS.Evec[2]
                  << std::endl;
        std::cout << "Test Z dot : " << test_RHS.Z << std::endl;
        std::cout << "My Z dot :   " << my_RHS.Z << std::endl;
    }

    return 0;
}
