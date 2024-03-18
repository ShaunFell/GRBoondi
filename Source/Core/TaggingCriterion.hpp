/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

 /*
 Code retrieved from 
 https://github.com/KAClough/GRChombo_public/blob/main/Examples/Oscillaton/HamTaggingCriterion.hpp

 Extraction piece from ChiExtractionTaggingCriterion.hpp
 */

#ifndef CUSTOMTAGGINGCRITERION_HPP_
#define CUSTOMTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"

class CustomTaggingCriterion
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


  public:
    CustomTaggingCriterion(double dx, const int a_level, const double a_L, const double a_rat,
                                  const std::array<double, CH_SPACEDIM> a_center,
                                  const spherical_extraction_params_t a_params,
                                  const bool a_activate_extraction = false,
                                  const bool a_activate_ham_tagging = false) : m_center{a_center}, m_dx(dx), m_L{a_L}, m_init_ratio{a_rat},m_params{a_params}, m_level{a_level}, m_activate_extraction{a_activate_extraction}, m_activate_ham_tagging{a_activate_ham_tagging}
                                  {
                                  };

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

    template <class data_t>
    data_t ExtractionTagging(Cell<data_t> current_cell) const
    {
        data_t criterion { 0. };

        if (m_activate_extraction)
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



    template <class data_t> 
    void compute(Cell<data_t> current_cell) const
    {
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);

        //first run fixed grid
        data_t FixedGridsTaggingCriterion { FixedGridTagging(current_cell) };

        //then run Extraction tagging
        data_t ExtractionTaggingCriterion { ExtractionTagging(current_cell) };

        data_t criterion { simd_max(FixedGridsTaggingCriterion, ExtractionTaggingCriterion) };

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* CUSTOMTAGGINGCRITERION_HPP_ */
