#ifndef CUSTOMTAGGINGCRITERION_HPP_
#define CUSTOMTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

class CustomTaggingCriterion
{
  protected:
    const double m_dx;
    const double m_L;
    const double m_init_ratio;
    const std::array<double, CH_SPACEDIM> m_center;
    const int m_level;
    const bool m_activate_gnn_tagging;
    const double m_threshold;

  public:
    CustomTaggingCriterion(double dx, int a_level, double a_L, double a_ratio,
                           double a_threshold,
                           std::array<double, CH_SPACEDIM> a_center,
                           bool activate_gnn_tagging)
        : m_dx(dx), m_L(a_L), m_init_ratio(a_ratio), m_center(a_center),
          m_level(a_level), m_activate_gnn_tagging(activate_gnn_tagging),
          m_threshold(a_threshold) {};

    template <class data_t>
    data_t FixedGridTagging(Cell<data_t> current_cell) const
    {
        data_t criterion{0.};
        // make sure the inner part is regridded around the horizon
        // m_init_ratio sets the initial ratio of the 1st refinement level
        // compared to the grid size m_init_ratio = 0.25 corresponds to tagging
        // cells around center +/ - L/4 Use the taxi cab metric to set the
        // regridding ratio
        double ratio = m_init_ratio * pow(2.0, -(m_level));
        const Coordinates<data_t> coords(current_cell, m_dx, m_center);
        const data_t max_abs_xy = simd_max(abs(coords.x), abs(coords.y));
        const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));
        auto regrid = simd_compare_lt(max_abs_xyz, m_L * ratio);
        criterion = simd_conditional(regrid, 100.0, criterion);

        return criterion;
    }

    template <class data_t> data_t gnnTagging(Cell<data_t> current_cell) const
    {
        // set default value
        data_t criterion{0.};

        if (m_activate_gnn_tagging)
        {
            auto gnn_data = current_cell.load_vars(
                c_gnn); // load the gnn diagnostic data from the grid
            auto gnn_data_abs =
                abs(gnn_data); // compute the absolute value, taking into
                               // account SIMD functionality
            auto gnn_criterion =
                m_threshold * m_threshold /
                gnn_data_abs; // Invert the gnn value to obtain the criterion
            criterion += gnn_criterion; // add back to the criterion value
        };

        return criterion;
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // compute each individual tagging criterion
        data_t FixedGridsTaggingCriterion{FixedGridTagging(current_cell)};
        data_t gnnTaggingCriterion{gnnTagging(current_cell)};

        // Compute the maximum of the two criteria
        data_t criterion{
            simd_max(FixedGridsTaggingCriterion, gnnTaggingCriterion)};

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif // CUSTOMTAGGINGCRITERION_HPP_