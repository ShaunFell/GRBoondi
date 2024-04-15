
// GRBoondi
#include "DataContainer.hpp"
#include "DataManipulation.hpp"
#include "SetupFunctions.hpp" //For setting up MPI processes
#include "SimpleDataReader.hpp"

#define SMALL_ERR 1e-4

int main(int argc, char *argv[])
{
    mainSetup(argc, argv); // setup MPI processes

    SimpleDataReader<double> reader{"KerrdeSitter_rPlus.dat"};
    DataContainer<double> data = reader.get_data();

    // Test the nearest neighbor algorithm
    std::cout << "###########################################" << std::endl;
    std::cout << "Testing nearest neighbor algorithm" << std::endl;
    std::cout << "###########################################" << std::endl;

    std::vector<std::vector<double>> coords{data.get_coords()};
    std::vector<double> vals{data.get_data()};

    std::vector<double> query_point{0.0526, 0.055};

    std::cout << std::setprecision(15);

    std::pair<std::vector<double>, std::vector<std::vector<double>>>
        nearest_3_neighbors{
            DataManipulation::find_nearest_neighbors(coords, query_point, 3)};
    std::cout << "Indicies of neighbors: " << std::endl;
    for (int i{0}; i < nearest_3_neighbors.first.size(); i++)
    {
        std::cout << nearest_3_neighbors.first[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Coords of neighbors: " << std::endl;
    for (int i{0}; i < nearest_3_neighbors.second.size(); i++)
    {
        for (int j{0}; j < nearest_3_neighbors.second[i].size(); j++)
        {
            std::cout << nearest_3_neighbors.second[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Test the interpolation algorithm
    std::cout << "###########################################" << std::endl;
    std::cout << "Testing interpolation algorithm" << std::endl;
    std::cout << "###########################################" << std::endl;

    // create 3D point
    std::vector<double> point1{nearest_3_neighbors.second[0][0],
                               nearest_3_neighbors.second[0][1],
                               vals[nearest_3_neighbors.first[0]]};
    std::vector<double> point2{nearest_3_neighbors.second[1][0],
                               nearest_3_neighbors.second[1][1],
                               vals[nearest_3_neighbors.first[1]]};
    std::vector<double> point3{nearest_3_neighbors.second[2][0],
                               nearest_3_neighbors.second[2][1],
                               vals[nearest_3_neighbors.first[2]]};

    std::cout << "Point 1: " << point1[0] << " " << point1[1] << " "
              << point1[2] << std::endl;
    std::cout << "Point 2: " << point2[0] << " " << point2[1] << " "
              << point2[2] << std::endl;
    std::cout << "Point 3: " << point3[0] << " " << point3[1] << " "
              << point3[2] << std::endl;

    std::cout << "Query Point: " << query_point[0] << " " << query_point[1]
              << std::endl;

    std::cout << "Interpolated value: "
              << DataManipulation::lin_interp_2d(point1, point2, point3,
                                                 query_point)
              << std::endl;

    mainFinalize(); // cleanup MPI processes
    return 0.;
}