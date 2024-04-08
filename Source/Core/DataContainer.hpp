/* GRBoondi
 * Please refer to LICENSE in GRBoondi's root directory.
 */
#ifndef DATACONTAINER_HPP_INCLUDED
#define DATACONTAINER_HPP_INCLUDED

/*
This class stores the calculated simulation data
This is a Value Container that stores the values explicitly and that
    doesn't allocate any new memory space
*/
template <typename data_type>
class DataContainer
{
    private:
        //vector of pairs, first member is the time data and second the derived data
        std::vector<std::pair<double, std::vector<data_type>>> m_data_history; 

        //current length of m_data_history
        int  m_current_size;

    public:

        DataContainer(): m_current_size{0}, m_data_history{} {}; //default constructor

        //constructor to initialize container with a single data point
        DataContainer(double time, std::vector<data_type> data): m_current_size{1}, m_data_history{{time, data}} {}; 

        //method to add a single new data point
        void update(double time, std::vector<data_type> data_point)
        {
            m_data_history.push_back({time, data_point});
            m_current_size += 1;
        }

        //method to grab array of time coordinates
        std::vector<double> get_time_data() const
        {
            std::vector<double> time_data;
            for (int i = 0; i < m_current_size; i++)
            {
                time_data.push_back(m_data_history[i].first);
            }
            return time_data;
        }

        //Overload operator[] to grab a particular data point from a given time coordinate
        std::vector<data_type>& operator[](double time)
        {
            //we iterate over the list of data points (in reverse) to find the time coordinate
            for (int i { m_current_size - 1 }; i >= 0; i--)
            {
                if (m_data_history[i].first == time)
                {
                    return m_data_history[i].second;
                }
            }
            //should never get here
            return m_data_history[0].second;
        }


};

#endif /* DATACONTAINER_HPP_INCLUDED */