/* GRBoondi
 * Please refer to LICENSE in GRBoondi's root directory.
 */
#ifndef DATACONTAINER_HPP_INCLUDED
#define DATACONTAINER_HPP_INCLUDED

#include <vector>

/*
This class stores the calculated simulation data
This is a Value Container that stores the values explicitly and that
    doesn't allocate any new memory space
*/
template <typename data_type> class TimeDataContainer
{
  private:
    // vector of pairs, first member is the time data and second the derived
    // data
    std::vector<std::pair<data_type, std::vector<data_type>>> m_data_history;

    // current length of m_data_history
    int m_current_size;

  public:
    TimeDataContainer()
        : m_current_size{0}, m_data_history{} {}; // default constructor

    // constructor to initialize container with a single data point
    TimeDataContainer(data_type time, std::vector<data_type> data)
        : m_current_size{1}, m_data_history{{time, data}} {};

    // method to add a single new data point
    void update(data_type time, std::vector<data_type> data_point)
    {
        m_data_history.push_back({time, data_point});
        m_current_size += 1;
    }

    // method to grab array of time coordinates
    std::vector<data_type> get_time_data() const
    {
        std::vector<data_type> time_data;
        for (int i = 0; i < m_current_size; i++)
        {
            time_data.push_back(m_data_history[i].first);
        }
        return time_data;
    }

    // Overload operator[] to grab a particular data point from a given time
    // coordinate
    std::vector<data_type> &operator[](data_type time)
    {
        // we iterate over the list of data points (in reverse) to find the time
        // coordinate
        for (int i{m_current_size - 1}; i >= 0; i--)
        {
            if (m_data_history[i].first == time)
            {
                return m_data_history[i].second;
            }
        }
        // should never get here
        return m_data_history[0].second;
    }
};

/*
This class stores generic N-dimensional scalar data, in the form of
coordinate-value pairs This is a value container that stores the values
explicitly and doesnt allocate new memory space
*/
template <typename data_type> class DataContainer
{
  private:
    // vector of pairs, first member is the coordinate location and the second
    // the scalar data
    std::vector<std::pair<std::vector<data_type>, data_type>> m_spatial_data;

    // current length of data
    int m_current_size;

  public:
    DataContainer()
        : m_current_size{0}, m_spatial_data{} {}; // default constructor

    // constructor to initialize container with an NxM matrix, assuming the last
    // column
    //  is the value data
    //  Sorting of input matrix is inherited
    DataContainer(std::vector<std::vector<data_type>> a_matrix_data);

    // Grab a particular data point by its index in the sorted data
    std::pair<std::vector<data_type>, data_type> &get_inx(int inx);

    // Get current size of data container
    int size() const;

    // add data to container
    void add_data(std::vector<data_type> coord, data_type data_point);

    // get nth column
    std::vector<data_type> get_column(int n);

    // get multiple columns
    std::vector<std::vector<data_type>> get_columns(int min_inx, int max_inx);

    // get matrix of coordinate data (the first n-1 columns)
    std::vector<std::vector<data_type>> get_coords();

    // get all the values in the container
    std::vector<data_type> get_data();
};

#include "DataContainer.impl.hpp"

#endif /* DATACONTAINER_HPP_INCLUDED */