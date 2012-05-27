//---------------------------------------------------------------------------//
/*!
 * \file Partitioner.hpp
 * \author Stuart R. Slattery
 * \brief Mesh partitioner unit tests.
 */
//---------------------------------------------------------------------------//

#include <algorithm>

#include "Partitioner.hpp"
#include "Mesh.hpp"

#include <UnitTestHelpers.hpp>

#include <Teuchos_ParameterList.hpp>

//---------------------------------------------------------------------------//
// Tests.
//---------------------------------------------------------------------------//

namespace Chimera
{

//---------------------------------------------------------------------------//
// Uniform mesh test.
TEUCHOS_UNIT_TEST( Partitioner, uniform_test )
{
    // Get the communicator.
    UnitTestHelpers::RCP_Comm comm = UnitTestHelpers::getDefaultComm();
    unsigned int my_rank = comm->getRank();
    unsigned int my_size = comm->getSize();

    // Block setup.
    int i_blocks = 1;
    int j_blocks = 1;
    if ( my_size == 2 )
    {
	j_blocks = 2;
    }
    else if ( my_size == 4 )
    {
	i_blocks = 2;
	j_blocks = 2;
    }
    double x_min = 0.0;
    double y_min = x_min;
    double x_max = 10.0*my_size;
    double y_max = x_max;
    int global_num_cells = 10*my_size;
    int x_local_num_cells = global_num_cells / i_blocks;
    int y_local_num_cells = global_num_cells / j_blocks;
    double cell_width = x_max / global_num_cells;

    // Setup a parameter list.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList );
    plist->set( "I_BLOCKS", i_blocks );
    plist->set( "J_BLOCKS", j_blocks );
    plist->set( "GRID_TYPE", "UNIFORM" );
    plist->set( "X_MIN", x_min );
    plist->set( "X_MAX", x_max );
    plist->set( "Y_MIN", y_min );
    plist->set( "Y_MAX", y_max );
    plist->set( "X_NUM_CELLS", global_num_cells );
    plist->set( "Y_NUM_CELLS", global_num_cells );

    // Partition.
    Partitioner partitioner( comm, plist );

    // Check the partitioning.
    if ( my_size == 1 )
    {
	TEST_ASSERT( partitioner.getNumBlocks().first == 1 );
	TEST_ASSERT( partitioner.getNumBlocks().second == 1 );
    }
    else if ( my_size == 2 )
    {
	TEST_ASSERT( partitioner.getNumBlocks().first == 1 );
	TEST_ASSERT( partitioner.getNumBlocks().second == 2 );
    }
    else if ( my_size == 4 )
    {
	TEST_ASSERT( partitioner.getNumBlocks().first == 2 );
	TEST_ASSERT( partitioner.getNumBlocks().second == 2 );
    }

    // Check the mesh.
    Partitioner::RCP_Mesh mesh = partitioner.getMesh();

    unsigned int my_j_block = 
	std::floor( my_rank / partitioner.getNumBlocks().first );
    unsigned int my_i_block = 
	my_rank - partitioner.getNumBlocks().first*my_j_block;

    TEST_ASSERT( (int) mesh->getGlobalNumCells().first == global_num_cells );
    TEST_ASSERT( (int) mesh->getGlobalNumCells().second == global_num_cells );
    TEST_ASSERT( mesh->getGlobalLowerBounds().first == x_min );
    TEST_ASSERT( mesh->getGlobalLowerBounds().second == y_min );
    TEST_ASSERT( mesh->getGlobalUpperBounds().first == x_max );
    TEST_ASSERT( mesh->getGlobalUpperBounds().second == y_max );
    TEST_ASSERT( mesh->getBlockID().first == my_i_block );
    TEST_ASSERT( mesh->getBlockID().second == my_j_block );
    TEST_ASSERT( (int) mesh->getLocalNumCells().first == x_local_num_cells );
    TEST_ASSERT( (int) mesh->getLocalNumCells().second == y_local_num_cells );

    std::vector<double>::const_iterator i_vec_it, j_vec_it;

    double x_edge_val = 0.0;
    for ( i_vec_it = mesh->getLocalEdges().first.begin();
	  i_vec_it != mesh->getLocalEdges().first.end();
	  ++i_vec_it )
    {
	TEST_ASSERT( *i_vec_it == 
		     my_rank*x_local_num_cells*cell_width + x_edge_val );
	x_edge_val += cell_width;
    }

    double y_edge_val = 0.0;
    for ( j_vec_it = mesh->getLocalEdges().second.begin();
	  j_vec_it != mesh->getLocalEdges().second.end();
	  ++j_vec_it )
    {
	TEST_ASSERT( *j_vec_it == 
		     my_rank*y_local_num_cells*cell_width + y_edge_val );
	y_edge_val += cell_width;
    }

    for ( i_vec_it = mesh->getLocalWidths().first.begin();
	  i_vec_it != mesh->getLocalWidths().first.end();
	  ++i_vec_it )
    {
	TEST_ASSERT( *i_vec_it == cell_width );
    }

    for ( j_vec_it = mesh->getLocalWidths().second.begin();
	  j_vec_it != mesh->getLocalWidths().second.end();
	  ++j_vec_it )
    {
	TEST_ASSERT( *j_vec_it == cell_width );
    }

    x_edge_val = 0.0;
    for ( i_vec_it = partitioner.getGlobalEdges().first.begin();
	  i_vec_it != partitioner.getGlobalEdges().first.end();
	  ++i_vec_it )
    {
	TEST_ASSERT( *i_vec_it == x_edge_val );
	x_edge_val += 1.0;
    }

    y_edge_val = 0.0;
    for ( j_vec_it = partitioner.getGlobalEdges().second.begin();
	  j_vec_it != partitioner.getGlobalEdges().second.end();
	  ++j_vec_it )
    {
	TEST_ASSERT( *j_vec_it == y_edge_val );
	y_edge_val += 1.0;
    }
}

//---------------------------------------------------------------------------//
// Nonuniform mesh test.
// TEUCHOS_UNIT_TEST( Partitioner, nonuniform_test )
// {
//     // Get the communicator.
//     UnitTestHelpers::RCP_Comm comm = UnitTestHelpers::getDefaultComm();
//     unsigned int my_rank = comm->getRank();
//     unsigned int my_size = comm->getSize();

//     // Block setup.
//     int i_blocks = 1;
//     int j_blocks = 1;
//     if ( my_size == 2 )
//     {
// 	j_blocks = 2;
//     }
//     else if ( my_size == 4 )
//     {
// 	i_blocks = 2;
// 	j_blocks = 2;
//     }
//     double x_min = 0.0;
//     double y_min = x_min;
//     double x_max = 10.0*my_size;
//     double y_max = x_max;
//     int global_num_cells = 10*my_size;
//     int x_local_num_cells = global_num_cells / i_blocks;
//     int y_local_num_cells = global_num_cells / j_blocks;
//     double cell_width = x_max / global_num_cells;

//     // Setup edge vectors.
//     std::vector<double> x_edges, y_edges;
//     double x_edge_val = 0.0;
//     for ( int i = 0; i < x_local_num_cells; ++i )
//     {
// 	x_edges.push_back( 
// 	    my_rank*cell_width*x_local_num_cells + x_edge_val );
// 	x_edge_val += cell_width;
//     }

//     double y_edge_val = 0.0;
//     for ( int j = 0; j < y_local_num_cells; ++j )
//     {
// 	y_edges.push_back( 
// 	    my_rank*cell_width*y_local_num_cells + y_edge_val );
// 	y_edge_val += cell_width;
//     }

//     // Setup a parameter list.
//     Teuchos::RCP<Teuchos::ParameterList> plist =
// 	Teuchos::rcp( new Teuchos::ParameterList );
//     plist->set( "I_BLOCKS", i_blocks );
//     plist->set( "J_BLOCKS", j_blocks );
//     plist->set( "GRID_TYPE", "NONUNIFORM" );
//     plist->set( "X_EDGES", x_edges );
//     plist->set( "Y_EDGES", y_edges );

//     // Partition.
//     Partitioner partitioner( comm, plist );

//     // Check the partitioning.
//     if ( my_size == 1 )
//     {
// 	TEST_ASSERT( partitioner.getNumBlocks().first == 1 );
// 	TEST_ASSERT( partitioner.getNumBlocks().second == 1 );
//     }
//     else if ( my_size == 2 )
//     {
// 	TEST_ASSERT( partitioner.getNumBlocks().first == 1 );
// 	TEST_ASSERT( partitioner.getNumBlocks().second == 2 );
//     }
//     else if ( my_size == 4 )
//     {
// 	TEST_ASSERT( partitioner.getNumBlocks().first == 2 );
// 	TEST_ASSERT( partitioner.getNumBlocks().second == 2 );
//     }

//     // Check the mesh.
//     Partitioner::RCP_Mesh mesh = partitioner.getMesh();

//     unsigned int my_j_block = 
// 	std::floor( my_rank / partitioner.getNumBlocks().first );
//     unsigned int my_i_block = 
// 	my_rank - partitioner.getNumBlocks().first*my_j_block;

//     TEST_ASSERT( (int) mesh->getGlobalNumCells().first == global_num_cells );
//     TEST_ASSERT( (int) mesh->getGlobalNumCells().second == global_num_cells );
//     TEST_ASSERT( mesh->getGlobalLowerBounds().first == x_min );
//     TEST_ASSERT( mesh->getGlobalLowerBounds().second == y_min );
//     TEST_ASSERT( mesh->getGlobalUpperBounds().first == x_max );
//     TEST_ASSERT( mesh->getGlobalUpperBounds().second == y_max );
//     TEST_ASSERT( mesh->getBlockID().first == my_i_block );
//     TEST_ASSERT( mesh->getBlockID().second == my_j_block );
//     TEST_ASSERT( (int) mesh->getLocalNumCells().first == x_local_num_cells );
//     TEST_ASSERT( (int) mesh->getLocalNumCells().second == y_local_num_cells );

//     std::vector<double>::const_iterator i_vec_it, j_vec_it;

//     x_edge_val = 0.0;
//     for ( i_vec_it = mesh->getLocalEdges().first.begin();
// 	  i_vec_it != mesh->getLocalEdges().first.end();
// 	  ++i_vec_it )
//     {
// 	TEST_ASSERT( *i_vec_it == 
// 		     my_rank*x_local_num_cells*cell_width + x_edge_val );
// 	x_edge_val += cell_width;
//     }

//     y_edge_val = 0.0;
//     for ( j_vec_it = mesh->getLocalEdges().first.begin();
// 	  j_vec_it != mesh->getLocalEdges().first.end();
// 	  ++j_vec_it )
//     {
// 	TEST_ASSERT( *j_vec_it == 
// 		     my_rank*y_local_num_cells*cell_width + y_edge_val );
// 	y_edge_val += cell_width;
//     }

//     for ( i_vec_it = mesh->getLocalWidths().first.begin();
// 	  i_vec_it != mesh->getLocalWidths().first.end();
// 	  ++i_vec_it )
//     {
// 	TEST_ASSERT( *i_vec_it == cell_width );
//     }

//     for ( j_vec_it = mesh->getLocalWidths().first.begin();
// 	  j_vec_it != mesh->getLocalWidths().first.end();
// 	  ++j_vec_it )
//     {
// 	TEST_ASSERT( *j_vec_it == cell_width );
//     }
// }

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Partitioner.hpp
//---------------------------------------------------------------------------//

