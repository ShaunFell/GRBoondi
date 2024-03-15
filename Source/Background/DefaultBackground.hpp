#ifndef DEFAULTBACKGROUND_HPP_
#define DEFAULTBACKGROUND_HPP_

#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"


class DefaultBackground
{
public:
    DefaultBackground::DefaultBackground() = default;

    template <class data_t>
    void compute(Cell<data_t> current_cell) const
    {
        //Set spacetime background to Minkowski in Cartesian coordinates

        vars.lapse = 1.;
        vars.K = 0.;

        FOR1(i)
        {
            vars.shift[i] = 0.;
            vars.d1_lapse[i] = 0.;

            FOR1(j)
            {
                vars.d1_shift[i][j] = 0.;
                vars.K_tensor[i][j] = 0.;
                if (i==j)
                {
                    vars.gamma[i][j] = 1.;
                } else 
                {
                    vars.gamma[i][j] = 0.;
                }
                FOR1(k)
                {
                    vars.d1_gamma[i][j][k] = 0.;
                }
            }
        }
    }

    bool check_if_excised(const Coordinates<double> &coords) const
    {
        return false; //Dont ever excise 
    }


}


#endif /* DEFAULTBACKGROUND_HPP_ */