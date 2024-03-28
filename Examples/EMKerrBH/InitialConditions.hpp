#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

/* 
Example class to calculate initial data
This example sets the initial data to be a uniform magnetic field
*/

#include "BaseProcaField.hpp"
#include "ADMProcaVars.hpp"
#include "KerrSchild.hpp" //background class
#include "Cell.hpp"
#include "TensorAlgebra.hpp"
#include "Tensor.hpp"


class Initial_EM_Conditions
{

    public:

        struct params_t
        {
            double init_amplitude;
        };

    protected:
        
        const params_t m_params;

        template <class data_t>
        using MatterVars = typename ADMProcaVars::MatterVars<data_t>;

    public:
        Initial_EM_Conditions(params_t a_params): m_params{a_params}{};

        template <class data_t>
        void compute(Cell<data_t> current_cell) const
        {
            MatterVars<data_t> matter_vars;;
            VarsTools::assign(matter_vars,0.); //assign all matter variables in this cell to zero;

            matter_vars.Avec[2] = m_params.init_amplitude; //assign z component of vector to amplitude;

            current_cell.store_vars(matter_vars); //push matter vars to cell
        }
};


#endif //INITIALCONDITIONS_HPP_