/*
GRBoondi 2024
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

#ifndef KERRSCHILDNEW_HPP_
#define KERRSCHILDNEW_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

#include "KerrSchild.hpp"


class KerrSchildNew : public KerrSchild {
  
  public:

    using KerrSchild::KerrSchild; //inherit constructor

    //overload the excision method to take in a 'buffer' parameter,
    // specifying the fractional distance from the horizon to excise
    bool check_if_excised(const Coordinates<double>& coords, const double buffer) const
    {
        // black hole params - mass M and spin a
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double a2 = a * a;

        // work out where we are on the grid
        const double x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        const double r_plus = M + sqrt(M * M - a2);
        const double r_minus = M - sqrt(M * M - a2);

        // position relative to outer horizon - 1 indicates on horizon
        // less than one is within
        const double outer_horizon =
            (x * x + y * y) / (2.0 * M * r_plus) + z * z / r_plus / r_plus;

        // position relative to inner horizon - 1 indicates on horizon, less
        // than 1 is within
        const double inner_horizon =
            (x * x + y * y) / (2.0 * M * r_minus) + z * z / r_minus / r_minus;

        bool is_excised = false;
        // value less than 1 indicates we are within the horizon
        if (outer_horizon < buffer || inner_horizon < 1/buffer)
        {
            is_excised = true;
        }
        return is_excised;
    }

};

#endif //KERRSCHILDNEW_HPP_ 