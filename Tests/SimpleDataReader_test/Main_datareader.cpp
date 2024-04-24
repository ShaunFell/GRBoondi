
// GRBoondi
#include "DataContainer.hpp"
#include "DataManipulation.hpp"
#include "SetupFunctions.hpp" //For setting up MPI processes
#include "SimpleDataReader.hpp"

#define SMALL_ERR 1e-4

// from
// https://stackoverflow.com/questions/17394149/how-can-you-efficiently-compare-two-stdvectors-for-equality-ignoring-the-orde
template <class T>
static bool compareVectors(const vector<T> &a, const vector<T> &b)
{
    const size_t n = a.size(); // make it const and unsigned!
    std::vector<bool> free(n, true);
    for (size_t i = 0; i < n; ++i)
    {
        bool matchFound = false;
        auto start = b.cbegin();
        while (true)
        {
            const auto position = std::find(start, b.cend(), a[i]);
            if (position == b.cend())
            {
                break; // nothing found
            }
            const auto index = position - b.cbegin();
            if (free[index])
            {
                // free pair found
                free[index] = false;
                matchFound = true;
                break;
            }
            else
            {
                start = position + 1; // search in the rest
            }
        }
        if (!matchFound)
        {
            return false;
        }
    }
    return true;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv); // setup MPI processes

    int reader_failed{0};
    int interp_failed{0};

    SimpleDataReader<double> reader{"data_test.dat_test"};
    DataContainer<double> data = reader.get_data();

    // Test the nearest neighbor algorithm
    std::cout << "###########################################" << std::endl;
    std::cout << "Testing nearest neighbor algorithm" << std::endl;
    std::cout << "###########################################" << std::endl;

    // get the coordinates and data values
    std::vector<std::vector<double>> coords{data.get_coords()};
    std::vector<double> vals{data.get_data()};

    // set the query point
    std::vector<double> query_point{2.5, 2.5};

    // find the nearest neighbors
    std::pair<std::vector<double>, std::vector<std::vector<double>>>
        nearest_neighbors{
            DataManipulation::find_nearest_neighbors(coords, query_point, 4)};

    // extract nearest neighbor coordinates
    std::vector<std::vector<double>> nearest_neighbor_coords{
        nearest_neighbors.second};

    // compare against known neighbors
    std::vector<std::vector<double>> ref_data{
        {2., 3.}, {3., 2.}, {3., 3.}, {2., 2.}};

    // compare the results
    if (compareVectors(ref_data, nearest_neighbor_coords))
    {
        std::cout << "Data reader test ..... PASSED" << std::endl;
    }
    else
    {
        std::cout << "Data reader test ..... FAILED" << std::endl;
        reader_failed = 1;
        std::cout << "nearest_neighbor_coords = ";
        for (auto el : nearest_neighbor_coords)
        {
            std::cout << "(";
            for (auto el2 : el)
            {
                std::cout << el2 << " ";
            }
            std::cout << ")";
        }
        std::cout << std::endl;
        std::cout << "ref_data = ";
        for (auto el : ref_data)
        {
            std::cout << "(";
            for (auto el2 : el)
            {
                std::cout << el2 << " ";
            }
            std::cout << ")";
        }
        std::cout << std::endl;
    }

    // Test the interpolation algorithm
    std::cout << "\n\n###########################################" << std::endl;
    std::cout << "Testing interpolation algorithm" << std::endl;
    std::cout << "###########################################" << std::endl;

    // create 3D point
    std::vector<double> point1{nearest_neighbors.second[0][0],
                               nearest_neighbors.second[0][1],
                               vals[nearest_neighbors.first[0]]};
    std::vector<double> point2{nearest_neighbors.second[1][0],
                               nearest_neighbors.second[1][1],
                               vals[nearest_neighbors.first[1]]};
    std::vector<double> point3{nearest_neighbors.second[2][0],
                               nearest_neighbors.second[2][1],
                               vals[nearest_neighbors.first[2]]};

    // get interpolated value
    double interp_data{
        DataManipulation::lin_interp_2d(point1, point2, point3, query_point)};

    // compare against known value
    double ref_interp_data{4.};
    if (interp_data == ref_interp_data)
    {
        std::cout << "Interpolation test ..... PASSED" << std::endl;
    }
    else
    {
        std::cout << "Interpolation test ..... FAILED" << std::endl;
        interp_failed = 1;
        std::cout << "interp_data = " << interp_data << std::endl;
        std::cout << "ref_interp_data = " << ref_interp_data << std::endl;
    }

    mainFinalize(); // cleanup MPI processes
    return (interp_failed || reader_failed);
}