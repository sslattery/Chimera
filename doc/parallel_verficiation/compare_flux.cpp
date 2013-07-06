#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <string>
#include <cassert>

int main( int argc, char * argv[] )
{
    std::string file_1 = "neutron.dat";
    std::string file_2 = "neutron.dat";

    int global_length;

    std::ifstream data_file_1;
    data_file_1.open( file_1.c_str() );
    data_file_1 >> global_length;
    std::vector<double> global_vector_1( global_length );
    for ( int i = 0; i < global_length; ++i )
    {
        data_file_1 >> global_vector_1[i];
    }
    data_file_1.close();

    std::ifstream data_file_2;
    data_file_2.open( file_2.c_str() );
    data_file_2 >> global_length;
    std::vector<double> global_vector_2( global_length );
    for ( int i = 0; i < global_length; ++i )
    {
        data_file_2 >> global_vector_2[i];
    }
    data_file_2.close();

    std::vector<double> difference( global_length );
    for ( int i = 0; i < global_length; ++i )
    {
        difference[i] = std::abs( global_vector_1[i] - global_vector_2[i] );
    }

    std::cout << "Max Difference: " 
              << *std::max_element( difference.begin(), difference.end() );

    std::cout << "Min Difference: " 
              << *std::min_element( difference.begin(), difference.end() );
}
