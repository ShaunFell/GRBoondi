#ifndef CUSTOMTAGGINGCRITERION_HPP_
#define CUSTOMTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "simd.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"


class CustomTaggingCriterion
{
    protected:
        const double m_dx;
        const double m_L;
        const double m_init_ratio;
        const std::array<double, CH_SPACEDIM> m_center;
        const int m_level;
        const bool m_activate_gnn_tagging;



    public:
        CustomTaggingCriterion(double dx, int a_level, double a_L, double a_ratio,
                                  std::array<double, CH_SPACEDIM> a_center, bool activate_gnn_tagging): m_dx(dx), m_L(a_L), m_init_ratio(a_ratio), m_center(a_center), m_level(a_level), m_activate_gnn_tagging(activate_gnn_tagging){};

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
    data_t gnnTagging(Cell<data_t> current_cell) const
    {
        data_t criterion { 0. };

        if (m_activate_gnn_tagging)
        {
            auto gnn_data = current_cell.load_vars(c_gnn);
            auto gnn_data_abs = abs(gnn_data);
            criterion += gnn_data_abs * m_dx;
        };

        return criterion;
    };

    template <class data_t>
    void compute(Cell<data_t> current_cell) const
    {
        data_t FixedGridsTaggingCriterion { FixedGridTagging(current_cell) };
        data_t gnnTaggingCriterion { gnnTagging(current_cell) };

        data_t criterion { simd_max(FixedGridsTaggingCriterion, gnnTaggingCriterion) };

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};



#endif // CUSTOMTAGGINGCRITERION_HPP_ 