
#ifndef ADMPROCAVARS_HPP_
#define ADMPROCAVARS_HPP_

#include "Tensor.hpp"
#include "VarsTools.hpp" //For mapping between Vars and Chombo grid

namespace ADMProcaVars
{
    template <class data_t> 
    struct MatterVars{
        data_t phi; //scalar part of Proca field
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
    }; //end of struct MatterVars

    //structure holding the matter field variables that require 2nd derivatives
    template <class data_t>
    struct Diff2MatterVars {
        Tensor<1, data_t> Avec;

        //provide function that maps between above Vars and Chombo grid variables
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function){
            using namespace VarsTools; //define_enum_mapping is part of VarsTools namespace
            define_enum_mapping(mapping_function, GRInterval<c_Avec1, c_Avec3>(), Avec);
        }
    }; //end of struct Diff2MatterVars
};

#endif //ADMPROCAVARS_HPP_