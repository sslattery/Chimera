//---------------------------------------------------------------------------//
/*!
 * \file tstMesh.cpp
 * \author Stuart R. Slattery
 * \brief Mesh unit tests.
 */
//---------------------------------------------------------------------------//

#include "Mesh.hpp"

#include <UnitTestHelpers.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//

namespace Chimera
{

TEUCHOS_UNIT_TEST( Mesh, parallel_mesh_test )
{
    // Get the communicator.
    UnitTestHelpers::RCP_Comm comm = UnitTestHelpers::getDefaultComm();
    unsigned int myRank = comm->getRank();
    unsigned int mySize = comm->getSize();

    // Setup local edge vectors.
    int local_edge_size = 11;
    int cell_width = 1.0;
    double lower_bounds = 0.0;
    double upper_bounds = 10.0*mySize;
    std::vector<double> x_edges, y_edges;
    double edge_val = 0.0;
    for ( int i = 0; i < local_edge_size; ++i )
    {
	x_edges.push_back( myRank*(local_edge_size-1) + edge_val );
	y_edges.push_back( edge_val );
	edge_val += cell_width;
    }

    // Build the mesh.
    unsigned int global_size = mySize*(local_edge_size-1) + 1;
    Mesh mesh( global_size, global_size,
	       lower_bounds, lower_bounds,
	       upper_bounds, upper_bounds,
	       myRank, myRank,
	       x_edges, y_edges );

    // Check the mesh.
    TEST_ASSERT( mesh.getNumGlobalCells().first == global_size );
    TEST_ASSERT( mesh.getNumGlobalCells().second == global_size );
    TEST_ASSERT( mesh.getLowerBounds().first == lower_bounds );
    TEST_ASSERT( mesh.getLowerBounds().second == lower_bounds );
    TEST_ASSERT( mesh.getBlockID().first == myRank );
    TEST_ASSERT( mesh.getBlockID().second == myRank );

    for ( int i = 0; i < local_edge_size; ++i )
    {
	TEST_ASSERT( mesh.getLocalEdges().first[i] == x_edges[i] );
	TEST_ASSERT( mesh.getLocalEdges().second[i] == y_edges[i] );
    }
    for ( int i = 0; i < local_edge_size-1; ++i )
    {
	TEST_ASSERT( mesh.getLocalWidths().first[i] == cell_width );
	TEST_ASSERT( mesh.getLocalWidths().second[i] == cell_width );
    }
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end tstMesh.cpp
//---------------------------------------------------------------------------//

