/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

#if !defined(DATACONTAINER_HPP_INCLUDED)
#error "This file should only be included through DataContainer.hpp"
#endif

#ifndef DATACONTAINER_IMPL_H_INCLUDED
#define DATACONTAINER_IMPL_H_INCLUDED

template <class data_type>
DataContainer<data_type>::DataContainer(std::vector<std::vector<data_type>> a_matrix_data)
{
    m_current_size = a_matrix_data.size();
    for (auto pt: a_matrix_data)
    {
        //assume value is the last element
        std::vector<data_type> coords{pt.begin(), pt.end() - 1};   
        data_type val{pt.back()};
        m_spatial_data.push_back({coords, val});
    }
}

template <class data_type>
std::pair<std::vector<data_type>, data_type>& DataContainer<data_type>::get_inx(int inx)
{
    return  m_spatial_data[inx];
}


template <class data_type>
int DataContainer<data_type>::size() const
{
    return m_current_size;
}

template <class data_type>
void DataContainer<data_type>::add_data(std::vector<data_type> coord, data_type data_point)
{
    m_spatial_data.push_back({coord, data_point});
    m_current_size += 1;
}

template <class data_type>
std::vector<data_type> DataContainer<data_type>::get_column(int n)
{
    std::vector<data_type> column;
    for (int i = 0; i < m_current_size; i++)
    {
        column.push_back(m_spatial_data[i].first[n]);
    }
    return column;
}

template <class data_type>
std::vector<std::vector<data_type>> DataContainer<data_type>::get_columns(int min_inx, int max_inx)
{
    std::vector<std::vector<data_type>> columns;
    for (int i = min_inx; i <= max_inx; i++)
    {
        columns.push_back(get_column(i));
    }
    return columns;
}

template <class data_type>
std::vector<std::vector<data_type>> DataContainer<data_type>::get_coords()
{
    std::vector<std::vector<data_type>> coords;
    for (auto pair: m_spatial_data)
    {
        coords.push_back(pair.first);
    }
    return coords;
}

template <class data_type>
std::vector<data_type> DataContainer<data_type>::get_data()
{
    std::vector<data_type> data;
    for (auto pair: m_spatial_data)
    {
        data.push_back(pair.second);
    }
    return data;
}



#endif /* DATACONTAINER_IMPL_H_INCLUDED */