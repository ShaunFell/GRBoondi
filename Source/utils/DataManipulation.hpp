/*
GRBoondi
Copyright 2024, Shaun Fell
Please refer to LICENSE in GRBoondi's root directory
*/


/*
This namespace contains various functions to manipulate data, such as nearest neighbor finders
*/

#ifndef DATAMANIPULATION_HPP_
#define DATAMANIPULATION_HPP_

// Chombo includes
#include "CH_assert.H"

#include <limits>
#include <cmath>

namespace DataManipulation
{

    //Brute force search for N nearest neighbors
    // Returns a pair object, first position is the indicies of the neighbors, second are the neighbors itself
    template <class data_t>
    std::pair< std::vector<data_t>, std::vector<std::vector<data_t>>> find_nearest_neighbors(std::vector<std::vector<data_t>> a_matrix_data, std::vector<data_t> a_query_point, int a_n_neighbors)
    {
        //container for the nearest neighbors
        std::vector<std::vector<data_t>> neighbors;
        std::vector<data_t> inx_neighbors_con;
        

        //create vector of indicies to keep track of original location of neighbor coords
        std::vector<data_t> data_indicies;
        for (int inx { 0 }; inx < a_matrix_data.size(); inx++)
        {
            data_indicies.push_back(inx);
        }

        //find the nearest neighbors
        for ( int inx_neighbors { 0 }; inx_neighbors < a_n_neighbors; inx_neighbors++)
        {
            //minimum distance for this neighbor
            data_t min_dist = std::numeric_limits<data_t>::max();

            //index of this neighbor
            int min_inx;

            //iterate over the rows in the data
            for (int coord_vec_inx { 0 }; coord_vec_inx < a_matrix_data.size(); coord_vec_inx++)
            {
                //calculate the distance between each point and the query point

                //current coordinate vector in data
                std::vector<data_t> coord_vec = a_matrix_data[coord_vec_inx];

                data_t dist = 0.;
                for ( int coord { 0 }; coord < coord_vec.size(); coord++)
                {
                    dist += (coord_vec[coord] - a_query_point[coord]) * (coord_vec[coord] - a_query_point[coord]);
                }
                if (std::sqrt(dist) < min_dist)
                {
                    min_dist = std::sqrt(dist);
                    min_inx = coord_vec_inx;
                }
            }

            neighbors.push_back(a_matrix_data[min_inx]);
            inx_neighbors_con.push_back(data_indicies[min_inx]);

            //drop this neighbor from the copy of the data
            a_matrix_data.erase(a_matrix_data.begin()+min_inx);
            data_indicies.erase(data_indicies.begin()+min_inx); // drop corresponding index

            //then repeat for the next nearest neighbor
        }

        return {inx_neighbors_con, neighbors} ;

    }

    // Linearily interpolate between 3 points in 2D space
    template <class data_t>
    data_t lin_interp_2d(std::vector<data_t> P1, std::vector<data_t> P2, std::vector<data_t> P3, std::vector<data_t> query_point)
    {
        // First 2 slots in the P1, P2, P3 vectors are the coordinates of the three points, Last slot is the value at those points
        // P1 = { x1, y1, value1 }, etc. 
        // query_points = { x, y };
        // return is the interpolated value

        // check all input vectors have same dimension
        CH_assert( P1.size() == P2.size() && P1.size() == P3.size() );

        // Find the normal vector of the unique plane intersecting the three points
        std::vector<data_t> P1P2;
        std::vector<data_t> P1P3;

        // Find spanning vectors of the three points
        for (int i = 0; i < P1.size(); i++)
        {
            P1P2.push_back(P2[i] - P1[i]);
            P1P3.push_back(P3[i] - P1[i]);
        }

        // Take curl to find the normal vector
        std::vector<data_t> n = {P1P2[1]*P1P3[2] - P1P2[2]*P1P3[1], P1P2[2]*P1P3[0] - P1P2[0]*P1P3[2], P1P2[0]*P1P3[1] - P1P2[1]*P1P3[0]};

        // if curl is zero, then points are co-linear, and therefore cannot be interpolated in 3d space
        CH_assert(n[0] != 0 && n[1] != 0 && n[2] != 0);

        // Compute the interpolated value at the given query point
        data_t x1 { query_point[0] };
        data_t x2 { query_point[1] };
        
        data_t x0 { P1[0] };    
        data_t y0 { P1[1] };
        data_t z0 { P1[2] };

        //the interpolated value
        data_t interpz { - (x1 - x0) * n[0] / n[2] - (x2 - y0) * n[1] / n[2] + z0 };

        return interpz;
    }

};



#endif /* DATAMANIPULATION_HPP_ */