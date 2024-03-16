#ifndef PROCAFIELD_H_INCLUDED
#define PROCAFIELD_H_INCLUDED

#include "ADMFixedBGVars.hpp" //For metric variables
#include "FourthOrderDerivatives.hpp" //For calculating derivatives
#include "Tensor.hpp" //For performing tensorial operations
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //For user-defined variables (e.g. see EMKerrBH)
#include "VarsTools.hpp" //For mapping between Vars and Chombo grid
#include "simd.hpp" //for SIMD operations
#include "DefaultBackground.hpp" //Minkowski background as default



template <class background_t>
class BaseProcaField
{
    protected: 

        const background_t m_background;

    public:
        //constructor, inputs are matter params
        BaseProcaField(background_t a_background): m_background{a_background}{};

        template <class data_t> 
        using MetricVars = ADMFixedBGVars::Vars<data_t>;

        template <class data_t> 
        struct Vars{
            data_t phi;
            data_t Z; //auxilliary damping scalar
            Tensor<1, data_t> Avec; //Spatial part of Proca field
            Tensor<1, data_t> Evec; //Electric part of Proca field strength tensor

            //provide function that maps between above Vars and Chombo grid variables
            template <typename mapping_function_t>
            void enum_mapping(mapping_function_t mapping_function){
                using namespace VarsTools; //define_enum_mapping is part of VarsTools namespace
                define_enum_mapping(mapping_function, c_phi, phi);
                define_enum_mapping(mapping_function, c_Z, Z);
                define_enum_mapping(mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
                define_enum_mapping(mapping_function, GRInterval<c_Evec1, c_Evec3>(), Evec);
            }
        }; //end of struct Vars

        //structure holding the matter field variables that require 2nd derivatives
        template <class data_t>
        struct Diff2Vars {
            Tensor<1, data_t> Avec;

            //provide function that maps between above Vars and Chombo grid variables
            template <typename mapping_function_t>
            void enum_mapping(mapping_function_t mapping_function){
                using namespace VarsTools; //define_enum_mapping is part of VarsTools namespace
                define_enum_mapping(mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
            }
        }; //end of struct Diff2Vars

        
        /* NOTE:
        Base GRChombo uses templated versions of the following functions. Since we employ virtual functions for the modifications, we can no longer use templates.
        So we specify the form of the following methods and use the declared types already
        
        Note: We overload each function. One with SIMD support, and one without
        Is there a better way of doing this without overloading and repeating the code?

        */

        //declare the types

    
        //with and without SIMD vectorization
        using SIMD =  simd<double>;
        using NSIMD =  double;

        //matter variables with and without SIMD vectorization
        using SIMD_vars_t = Vars<SIMD>;
        using NSIMD_vars_t = Vars<NSIMD>;

        //derivatives of matter variables with and without SIMD vectorization
        using SIMD_vars_d1_t = Vars<Tensor<1, SIMD>>;
        using NSIMD_vars_d1_t = Vars<Tensor<1, NSIMD>>;

        //rhs variables with and without SIMD vectorization
        using SIMD_rhs_vars_t = Vars<SIMD>;
        using NSIMD_rhs_vars_t = Vars<NSIMD>;

        //metric variables with and without SIMD vectorization
        using SIMD_metric_vars_t = MetricVars<SIMD>;
        using NSIMD_metric_vars_t = MetricVars<NSIMD>;

        //second derivative of matter variables with and without SIMD vectorization
        using SIMD_diff2_vars_t = Diff2Vars<Tensor<2,SIMD>>;
        using NSIMD_diff2_vars_t = Diff2Vars<Tensor<2,NSIMD>>;


        //we put the method definitions in this header file to make use of above type aliases

        //with SIMD vectorization
        emtensor_t<SIMD> compute_emtensor(
            const SIMD_vars_t &matter_vars, //the value of the variables
            const SIMD_metric_vars_t &metric_vars,
            const SIMD_vars_d1_t &d1, //the 1st derivatives
            const Tensor<2, SIMD> &gamma_UU, //the inverse spatial metric
            const Tensor<3, SIMD> &chris_ULL //physical christoffel symbols
        ) const 
        {
            emtensor_t<SIMD> out;

            // D_i A_j  3-covariant derivative of spatial covector
            Tensor<2, SIMD> DA;
            FOR2(i, j)
            {
                DA[i][j] = d1.Avec[j][i];
                FOR1(k) { DA[i][j] -= chris_ULL[k][i][j] * matter_vars.Avec[k]; };
            };

            // D_i A_j - D_j A_i
            Tensor<2, SIMD> DA_antisym;
            FOR2(i, j) { DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; };

            // Electric Field Norm
            SIMD Enorm{0};
            FOR2(i, j) { Enorm += metric_vars.gamma[i][j] * matter_vars.Evec[i] * matter_vars.Evec[j]; };

            /////Components of EM tensor

            // Eulerian Energy //Checked. Agrees with Mathematica notebook
            out.rho = 1. / 2. * Enorm;
            FOR4(i, j, k, l)
            {
                out.rho += 1. / 2. * gamma_UU[k][i] * gamma_UU[l][j] * DA[i][j] *
                        DA_antisym[k][l];
            };

            // Eulerian Momentum //Checked. Agrees with Mathematica notebook.
            FOR1(i)
            {
                out.Si[i] = 0; // zero initialize
                FOR1(j) { out.Si[i] += matter_vars.Evec[j] * DA_antisym[i][j]; };
            };

            // Eulerian Stress //Chedked. Agrees with Mathematica notebook.
            FOR2(i, j)
            {
                out.Sij[i][j] = 0; // zero initialize

                out.Sij[i][j] += 1. / 2. * metric_vars.gamma[i][j] * Enorm;

                FOR2(l, k)
                {
                    out.Sij[i][j] +=
                        -metric_vars.gamma[i][l] * metric_vars.gamma[j][k] * matter_vars.Evec[l] * matter_vars.Evec[k] +
                        gamma_UU[k][l] * DA_antisym[i][l] * DA_antisym[j][k];

                    FOR2(m, n)
                    {
                        out.Sij[i][j] += -metric_vars.gamma[i][j] * gamma_UU[m][l] *
                                        gamma_UU[n][k] * DA[m][n] * DA_antisym[l][k];
                    };
                };
            };

            // Eulerian Stress scalar
            out.S = 0.0;
            FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; };

            //add modifications
            compute_emtensor_modification(out, matter_vars, metric_vars, d1, gamma_UU, chris_ULL);

            return out;
        };

        //without SIMD vectorization
        emtensor_t<NSIMD> compute_emtensor(
            const NSIMD_vars_t &matter_vars, //the value of the variables
            const NSIMD_metric_vars_t &metric_vars,
            const NSIMD_vars_d1_t &d1, //the 1st derivatives
            const Tensor<2, NSIMD> &gamma_UU, //the inverse spatial metric
            const Tensor<3, NSIMD> &chris_ULL //physical christoffel symbols
        ) const
        {

            emtensor_t<NSIMD> out;

            // D_i A_j  3-covariant derivative of spatial covector
            Tensor<2, NSIMD> DA;
            FOR2(i, j)
            {
                DA[i][j] = d1.Avec[j][i];
                FOR1(k) { DA[i][j] -= chris_ULL[k][i][j] * matter_vars.Avec[k]; };
            };

            // D_i A_j - D_j A_i
            Tensor<2, NSIMD> DA_antisym;
            FOR2(i, j) { DA_antisym[i][j] = d1.Avec[j][i] - d1.Avec[i][j]; };

            // Electric Field Norm
            NSIMD Enorm{0};
            FOR2(i, j) { Enorm += metric_vars.gamma[i][j] * matter_vars.Evec[i] * matter_vars.Evec[j]; };

            /////Components of EM tensor

            // Eulerian Energy //Checked. Agrees with Mathematica notebook
            out.rho = 1. / 2. * Enorm;
            FOR4(i, j, k, l)
            {
                out.rho += 1. / 2. * gamma_UU[k][i] * gamma_UU[l][j] * DA[i][j] *
                        DA_antisym[k][l];
            };

            // Eulerian Momentum //Checked. Agrees with Mathematica notebook.
            FOR1(i)
            {
                out.Si[i] = 0; // zero initialize
                FOR1(j) { out.Si[i] += matter_vars.Evec[j] * DA_antisym[i][j]; };
            };

            // Eulerian Stress //Chedked. Agrees with Mathematica notebook.
            FOR2(i, j)
            {
                out.Sij[i][j] = 0; // zero initialize

                out.Sij[i][j] += 1. / 2. * metric_vars.gamma[i][j] * Enorm;

                FOR2(l, k)
                {
                    out.Sij[i][j] +=
                        -metric_vars.gamma[i][l] * metric_vars.gamma[j][k] * matter_vars.Evec[l] * matter_vars.Evec[k] +
                        gamma_UU[k][l] * DA_antisym[i][l] * DA_antisym[j][k];

                    FOR2(m, n)
                    {
                        out.Sij[i][j] += -metric_vars.gamma[i][j] * gamma_UU[m][l] *
                                        gamma_UU[n][k] * DA[m][n] * DA_antisym[l][k];
                    };
                };
            };

            // Eulerian Stress scalar
            out.S = 0.0;
            FOR2(i, j) { out.S += out.Sij[i][j] * gamma_UU[i][j]; };

            //add modifications
            compute_emtensor_modification(out, matter_vars, metric_vars, d1, gamma_UU, chris_ULL);

            return out;
        };


        //with SIMD vectorization
        void matter_rhs(
            SIMD_rhs_vars_t &total_rhs, //RHS terms for all vars
            const SIMD_vars_t &matter_vars, //the value fo the variables
            const SIMD_metric_vars_t &metric_vars,
            const SIMD_vars_d1_t &d1, //the 1st derivatives
            const SIMD_diff2_vars_t &d2, //the 2nd derivatives
            const SIMD_vars_t &advec //value of the beta^i d_i(var) terms
        ) const
        {

            // calculate conformal contravariant metric and conformal christoffel
            // symbols
            const Tensor<2, SIMD> gamma_UU = TensorAlgebra::compute_inverse(metric_vars.gamma);
            const Tensor<3, SIMD> chris_phys =
                TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL;
            
            // covariant derivative of spatial part of Proca field
            Tensor<2, SIMD> DA;
            FOR2(i, j)
            {
                DA[i][j] = d1.Avec[j][i];
                FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
            }


            // evolution equations for spatial part of vector field (index down). 
            FOR1(i)
            {
                total_rhs.Avec[i] =
                    -metric_vars.lapse * d1.phi[i] - matter_vars.phi * metric_vars.d1_lapse[i] + advec.Avec[i];

                FOR1(j)
                {
                    total_rhs.Avec[i] += -metric_vars.lapse * metric_vars.gamma[i][j] * matter_vars.Evec[j] +
                                        matter_vars.Avec[j] * metric_vars.d1_shift[j][i];
                };
            };

            // evolution equations for Electric vector field (index up)
            FOR1(i)
            {
                total_rhs.Evec[i] = metric_vars.lapse * metric_vars.K * matter_vars.Evec[i] + advec.Evec[i];

                FOR1(j)
                {
                    total_rhs.Evec[i] += - matter_vars.Evec[j] * metric_vars.d1_shift[i][j];
                }

                FOR3(j, k, l)
                {
                    total_rhs.Evec[i] += gamma_UU[j][k] * gamma_UU[i][l] * (
                                            metric_vars.d1_lapse[j] * ( d1.Avec[k][l] - d1.Avec[l][k] ) +
                                            metric_vars.lapse * ( d2.Avec[k][l][j] - d2.Avec[l][k][j] )
                                        );

                    FOR1(m)
                    {
                        total_rhs.Evec[i] +=
                            -metric_vars.lapse * gamma_UU[j][k] * gamma_UU[i][l] *
                            (chris_phys[m][j][l] * (d1.Avec[k][m] - d1.Avec[m][k]) +
                            chris_phys[m][j][k] * (d1.Avec[m][l] - d1.Avec[l][m]));
                    };
                };
            };

            // evolution equation for auxiliary constraint-damping scalar field Z
            total_rhs.Z = 0.;

            //Evolution equation for phi could be determined from constraint equation
            total_rhs.phi = 0;

            //add modifications
            matter_rhs_modification(total_rhs, matter_vars, metric_vars, d1, d2, advec);
        }

        //without SIMD vectorization
        void matter_rhs(
            NSIMD_rhs_vars_t &total_rhs, //RHS terms for all vars
            const NSIMD_vars_t &matter_vars, //the value fo the variables
            const NSIMD_metric_vars_t &metric_vars,
            const NSIMD_vars_d1_t &d1, //the 1st derivatives
            const NSIMD_diff2_vars_t &d2, //the 2nd derivatives
            const NSIMD_vars_t &advec //value of the beta^i d_i(var) terms
        ) const
        {

            // calculate conformal contravariant metric and conformal christoffel
            // symbols
            const Tensor<2, NSIMD> gamma_UU = TensorAlgebra::compute_inverse(metric_vars.gamma);
            const Tensor<3, NSIMD> chris_phys =
                TensorAlgebra::compute_christoffel(metric_vars.d1_gamma, gamma_UU).ULL;
            
            // covariant derivative of spatial part of Proca field
            Tensor<2, NSIMD> DA;
            FOR2(i, j)
            {
                DA[i][j] = d1.Avec[j][i];
                FOR1(k) { DA[i][j] -= chris_phys[k][i][j] * matter_vars.Avec[k]; }
            }


            // evolution equations for spatial part of vector field (index down). 
            FOR1(i)
            {
                total_rhs.Avec[i] =
                    -metric_vars.lapse * d1.phi[i] - matter_vars.phi * metric_vars.d1_lapse[i] + advec.Avec[i];

                FOR1(j)
                {
                    total_rhs.Avec[i] += -metric_vars.lapse * metric_vars.gamma[i][j] * matter_vars.Evec[j] +
                                        matter_vars.Avec[j] * metric_vars.d1_shift[j][i];
                };
            };

            // evolution equations for Electric vector field (index up)
            FOR1(i)
            {
                total_rhs.Evec[i] = metric_vars.lapse * metric_vars.K * matter_vars.Evec[i] + advec.Evec[i];

                FOR1(j)
                {
                    total_rhs.Evec[i] += - matter_vars.Evec[j] * metric_vars.d1_shift[i][j];
                }

                FOR3(j, k, l)
                {
                    total_rhs.Evec[i] += gamma_UU[j][k] * gamma_UU[i][l] * (
                                            metric_vars.d1_lapse[j] * ( d1.Avec[k][l] - d1.Avec[l][k] ) +
                                            metric_vars.lapse * ( d2.Avec[k][l][j] - d2.Avec[l][k][j] )
                                        );

                    FOR1(m)
                    {
                        total_rhs.Evec[i] +=
                            -metric_vars.lapse * gamma_UU[j][k] * gamma_UU[i][l] *
                            (chris_phys[m][j][l] * (d1.Avec[k][m] - d1.Avec[m][k]) +
                            chris_phys[m][j][k] * (d1.Avec[m][l] - d1.Avec[l][m]));
                    };
                };
            };

            // evolution equation for auxiliary constraint-damping scalar field Z
            total_rhs.Z = 0.;

            //Evolution equation for phi could be determined from constraint equation
            total_rhs.phi = 0;

            //add modifications
            matter_rhs_modification(total_rhs, matter_vars, metric_vars, d1, d2, advec);
        }

        //With SIMD vectorization
        virtual void compute_emtensor_modification(
            emtensor_t<SIMD> &base_emtensor, //pass by reference to allow modifications
            const SIMD_vars_t &matter_vars,
            const SIMD_metric_vars_t &metric_vars,
            const SIMD_vars_d1_t &d1,
            const Tensor<2, SIMD> &gamma_UU,
            const Tensor<3, SIMD> &chris_ULL
        ) const {}; 

        //Without SIMD vectorization
        virtual void compute_emtensor_modification(
            emtensor_t<NSIMD> &base_emtensor,
            const NSIMD_vars_t &matter_vars,
            const NSIMD_metric_vars_t &metric_vars,
            const NSIMD_vars_d1_t &d1,
            const Tensor<2, NSIMD> &gamma_UU,
            const Tensor<3, NSIMD> &chris_ULL
        ) const {};

        //With SIMD vectorization
        virtual void matter_rhs_modification(
            SIMD_rhs_vars_t &total_rhs, //pass RHS by reference to allow modifications
            const SIMD_vars_t &matter_vars,
            const SIMD_metric_vars_t &metric_vars,
            const SIMD_vars_d1_t &d1,
            const SIMD_diff2_vars_t &d2,
            const SIMD_vars_t &advec
        ) const {};

        //Without SIMD vectorization
        virtual void matter_rhs_modification(
            NSIMD_rhs_vars_t &total_rhs,
            const NSIMD_vars_t &matter_vars,
            const NSIMD_metric_vars_t &metric_vars,
            const NSIMD_vars_d1_t &d1,
            const NSIMD_diff2_vars_t &d2,
            const NSIMD_vars_t &advec
        ) const {};

    
        //the templated version
        /*//method that computes EM tensor, given vars and derivatives
        template <class data_t, template <typename> class vars_t>
        emtensor_t<data_t> compute_emtensor(
            const vars_t<data_t> &matter_vars, //the value of the variables
            const MetricVars<data_t> &metric_vars,
            const vars_t<Tensor<1,data_t>> &d1, //the 1st derivatives
            const Tensor<2, data_t> &gamma_UU, //the inverse spatial metric
            const Tensor<3, data_t> &chris_ULL //physical christoffel symbols
        ) const;  */


        // the templated version
        //method which adds in the matter field RHS, given vars and derivatives
    /*  template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
        void matter_rhs(
            rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
            const vars_t<data_t> &vars, //the value fo the variables
            const MetricVars<data_t> &metric_vars,
            const vars_t<Tensor<1, data_t>> &d1, //value of 1st derivs
            const diff2_vars_t<Tensor<2, data_t>> &d2, //2nd derivs
            const vars_t<data_t> &advec //value of the beta^i d_i(var) terms
        ) const;  */

        //Virtual functions to allow modification of the matter field
    /*  template <class data_t, template <typename> class vars_t>
        virtual void compute_emtensor_modification(
            emtensor_t<data_t> base_emtensor,
            const vars_t<data_t> &matter_vars, //the value of the variables
            const MetricVars<data_t> &metric_vars,
            const vars_t<Tensor<1,data_t>> &d1, //the 1st derivatives
            const Tensor<2, data_t> &gamma_UU, //the inverse spatial metric
            const Tensor<3, data_t> &chris_ULL //physical christoffel symbols
        ) const {};  */
        
        /* template <class data_t, template <typename> class vars_t, template <typename> class diff2_vars_t, template <typename> class rhs_vars_t>
        virtual void matter_rhs_modification(
            rhs_vars_t<data_t> &total_rhs, //RHS terms for all vars
            const vars_t<data_t> &vars, //the value fo the variables
            const MetricVars<data_t> &metric_vars,
            const vars_t<Tensor<1, data_t>> &d1, //value of 1st derivs
            const diff2_vars_t<Tensor<2, data_t>> &d2, //2nd derivs
            const vars_t<data_t> &advec //value of the beta^i d_i(var) terms
        ) const {};  */

};


#include "BaseProcaField.impl.hpp"
#endif //PROCAFIELD_H_INCLUDED

