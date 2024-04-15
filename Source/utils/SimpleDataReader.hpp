/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/

#ifndef SIMPLEDATAREADER_HPP_
#define SIMPLEDATAREADER_HPP_

#include "DataContainer.hpp"

//GRChombo includes
#include "SmallDataIOReader.hpp"

// system includes
#include <fstream>
#include <string>
#include <vector>


// A class that reads in a simple data file, containing columns of data
// Inherits from SmallDataIOReader
// Adds ability to read in data to DataContainer
template <class data_type>
class SimpleDataReader: public SmallDataIOReader
{
    public:
        // Read in data during constructor call
        SimpleDataReader(std::string a_file_name)
        {
            open(a_file_name); // Open the file object
            determine_file_structure(); // Parse the file
        }

        ~SimpleDataReader() 
        {
            if (m_file.is_open()) 
            {
                m_file.close(); // close the file object if its open
            }
        }

        DataContainer<data_type> get_data()
        {
            // Parse all data in file
            std::vector<std::vector<data_type>> data { get_all_columns() };

            //transpose the matrix of column vectors
            std::vector<std::vector<data_type>>  transposed_data { transpose_data(data) };

            // Create data container object, and add data to object
            DataContainer<data_type> new_container{transposed_data};

            // Return data container
            return new_container;
        }


        std::vector<std::vector<data_type>> transpose_data(std::vector<std::vector<data_type>> data)
        {
            // Tranpose the vector of columns into a vector of rows
            std::vector<std::vector<data_type>> transposed_data;
            for (int i = 0; i < data[0].size(); i++)
            {
                std::vector<data_type> row;
                for (int j = 0; j < data.size(); j++)
                {
                    row.push_back(data[j][i]);
                }
                transposed_data.push_back(row);
            }
            return transposed_data;
        }

};







#endif /* SIMPLEDATAREADER_HPP_ */