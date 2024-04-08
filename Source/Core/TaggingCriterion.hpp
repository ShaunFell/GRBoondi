/* GRBoondi 2024
 * Please refer to LICENSE in GRBoondi's root directory.
 */

/*
This tagging criterion is a generic criterion that uses a fixed grid, but can also deal with flux extraction in the outer region of the simulation
*/

#ifndef TAGGINGCRITERION_HPP_
#define TAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

class TaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const double m_init_ratio;
    const spherical_extraction_params_t m_params;
    const std::array<double, CH_SPACEDIM> m_center;
    const int m_level;
    const bool m_activate_extraction;
    const bool m_activate_ham_tagging;
    const bool m_activate_extraction_tagging;
    const bool m_DIAGNOSTIC;


  public:
    TaggingCriterion(double dx, const int a_level, const double a_L, const double a_rat,
                                  const std::array<double, CH_SPACEDIM> a_center,
                                  const spherical_extraction_params_t a_params,
                                  const bool a_activate_extraction = false,
                                  const bool a_activate_ham_tagging = false,
                                  const bool a_activate_extraction_tagging = false,
                                  const bool DIAGNOSTIC = false) : m_center{a_center}, m_dx(dx), m_L{a_L}, m_init_ratio{a_rat},m_params{a_params}, m_level{a_level}, m_activate_extraction{a_activate_extraction}, m_activate_ham_tagging{a_activate_ham_tagging}, m_activate_extraction_tagging{a_activate_extraction_tagging}, m_DIAGNOSTIC{DIAGNOSTIC}
                                  {
                                  };

    //Fixed grid criteria. Grid levels are determined by their distance from the user-defined center of the grid.
    //Overall scaling and initial scaling are determined by the user
    template <class data_t> 
    data_t FixedGridTagging(Cell<data_t> current_cell) const
    {
        data_t criterion { 0. };
        // make sure the inner part is regridded around the horizon
        // take L as the length of full grid, so tag inner 1/2
        // of it, which means inner \pm L/4
        double ratio = m_init_ratio * pow(2.0, -(m_level));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
        auto regrid = simd_compare_lt(max_abs_xyz, m_L * ratio);
        criterion = simd_conditional(regrid, 100.0, criterion);

        return criterion;
    }

    //Extraction tagging. We want to have a fine enough resolution at the extraction surface
    template <class data_t>
    data_t ExtractionTagging(Cell<data_t> current_cell) const
    {
        data_t criterion { 0. };

        if (m_activate_extraction && m_activate_extraction_tagging)
        {
            const Coordinates<data_t> coords(current_cell, m_dx, m_center);
            
            for (int iradius {0}; iradius < m_params.num_extraction_radii; ++iradius)
            {
                //regrid if within extraction level
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const data_t r = coords.get_radius();
                    // add 20% buffer to extraction zone to avoid boundary
                    auto regrid = simd_compare_lt(r, 1.2*m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        return criterion;
    }


    //The compute function called by BoxLoops
    //Computes each criteria, then combines them into a single criterion and stores to grid
    template <class data_t> 
    void compute(Cell<data_t> current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);

        //first run fixed grid
        data_t FixedGridsTaggingCriterion { FixedGridTagging(current_cell) };

        //then run Extraction tagging
        data_t ExtractionTaggingCriterion { ExtractionTagging(current_cell) };

        //Combine criteria
        data_t criterion { simd_max(FixedGridsTaggingCriterion, ExtractionTaggingCriterion) };

        // Write back into the flattened Chombo box. If Diagnostic is turned on, then write to c_Tagging_Diagnostic

        current_cell.store_vars(criterion, 0);
  
    }
};

#endif /* TAGGINGCRITERION_HPP_ */
