//---------------------------------------------------------------------------//
/*!
 * \file Partitioner.cpp
 * \author Stuart R. Slattery
 * \brief Mesh partitioner defintion.
 */
//---------------------------------------------------------------------------//

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

#include "Partitioner.hpp"

#include <Chimera_Assertion.hpp>

#include <Teuchos_ENull.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Partitioner::Partitioner( const RCP_Comm &comm, const RCP_ParameterList &plist )
    : d_num_blocks( SizePair( plist->get<int>( "I_BLOCKS" ), 
			      plist->get<int>( "J_BLOCKS" ) ) )
{
    // Comm parameters.
    unsigned int my_rank = comm->getRank();
    unsigned int my_size = comm->getSize();

    // Check that the block specification and communicator are consistent.
    testAssertion( d_num_blocks.first * d_num_blocks.second == my_size );

    // Block indices.
    int my_j_block = std::floor( my_rank / d_num_blocks.first );
    int my_i_block = my_rank - d_num_blocks.first*my_j_block;

    std::vector<double> i_edges, j_edges;
    double global_i_min, global_i_max, global_j_min, global_j_max;
    int global_num_i, global_num_j;

    // Uniform grid case.
    if ( plist->get<std::string>("GRID_TYPE") == "UNIFORM" )
    {
	// Get the parameters.
	global_i_min = plist->get<double>( "X_MIN" );
	global_i_max = plist->get<double>( "X_MAX" );
	global_j_min = plist->get<double>( "Y_MIN" );
	global_j_max = plist->get<double>( "Y_MAX" );
	global_num_i = plist->get<int>( "X_NUM_CELLS" );
	global_num_j = plist->get<int>( "Y_NUM_CELLS" );

	// Cell widths.
	double width_i = (global_i_max - global_i_min) / global_num_i;
	double width_j = (global_j_max - global_j_min) / global_num_j;
	d_cell_size.first = width_i;
	d_cell_size.second = width_j;

	// Number of local vertices.
	int i_edges_size = std::floor( global_num_i / d_num_blocks.first ) + 1;
	int j_edges_size = std::floor( global_num_j / d_num_blocks.second ) + 1;

	// Remaining vertices to be tacked onto last blocks.
	int i_remainder = global_num_i % d_num_blocks.first;
	int j_remainder = global_num_j % d_num_blocks.second;

	// Set the local I edges.
	double i_edge_val = 0.0;
	for ( int i = 0; i < i_edges_size; ++i )
	{
	    i_edges.push_back( 
		my_i_block*width_i*(i_edges_size-1) + i_edge_val + global_i_min );

	    i_edge_val += width_i;
	}
	if ( my_i_block == (int) d_num_blocks.first )
	{
	    for ( int i = 0; i < i_remainder; ++i )
	    {
		i_edges.push_back( 
		    my_i_block*width_i*(i_edges_size-1) + i_edge_val 
		    + global_i_min );

		i_edge_val += width_i;
	    }
	}

	// Set the local J edges.
	double j_edge_val = 0.0;
	for ( int j = 0; j < j_edges_size; ++j )
	{
	    j_edges.push_back( 
		my_j_block*width_j*(j_edges_size-1) + j_edge_val + global_j_min );

	    j_edge_val += width_j;
	}
	if ( my_j_block == (int) d_num_blocks.second )
	{
	    for ( int j = 0; j < j_remainder; ++j )
	    {
		j_edges.push_back( 
		    my_j_block*width_j*(j_edges_size-1) + j_edge_val 
		    + global_j_min );

		j_edge_val += width_j;
	    }
	}

	// Set the global I edges.
	i_edge_val = 0.0;
	for ( int i = 0; i < global_num_i+1; ++i )
	{
	    d_global_edges.first.push_back( i_edge_val + global_i_min );
	    i_edge_val += width_i;
	}

	// Set the global I edges.
	j_edge_val = 0.0;
	for ( int j = 0; j < global_num_j+1; ++j )
	{
	    d_global_edges.second.push_back( j_edge_val + global_j_min );
	    j_edge_val += width_j;
	}
    }

    // Unsupported cases.
    else
    {
	testPrecondition( plist->get<std::string>("GRID_TYPE") == "UNIFORM" );
    }

    // Create th mesh.
    d_mesh = Teuchos::rcp( new Mesh( global_num_i, global_num_j,
				     global_i_min, global_j_min,
				     global_i_max, global_j_max,
				     my_i_block, my_j_block,
				     i_edges, j_edges ) );

    testPostcondition( d_mesh != Teuchos::null );

    // Set the local rows.
    int row_idx, idx_i, idx_j;
    for ( int j = 0; j < (int) j_edges.size()-1; ++j )
    {
	for ( int i = 0; i < (int) i_edges.size()-1; ++i )
	{
	    idx_i = my_i_block*d_mesh->getLocalNumCells().first + i;
	    idx_j = my_j_block*d_mesh->getLocalNumCells().second + j;
	    row_idx = idx_i + idx_j*d_global_edges.first.size();
	    d_local_rows.push_back( row_idx );
	}
    }
    if ( my_j_block == (int) d_num_blocks.second-1 )
    {
	int j = j_edges.size() - 1;

	for ( int i = 0; i < (int) i_edges.size()-1; ++i )
	{
	    idx_i = my_i_block*d_mesh->getLocalNumCells().first + i;
	    idx_j = my_j_block*d_mesh->getLocalNumCells().second + j;
	    row_idx = idx_i + idx_j*d_global_edges.first.size();
	    d_local_rows.push_back( row_idx );
	}
    }
    if ( my_i_block == (int) d_num_blocks.first-1 )
    {
	int i = i_edges.size() - 1;

	for ( int j = 0; j < (int) j_edges.size()-1; ++j )
	{
	    idx_i = my_i_block*d_mesh->getLocalNumCells().first + i;
	    idx_j = my_j_block*d_mesh->getLocalNumCells().second + j;
	    row_idx = idx_i + idx_j*d_global_edges.first.size();
	    d_local_rows.push_back( row_idx );
	}
    }
    if ( my_i_block == (int) d_num_blocks.first-1 &&
	 my_j_block == (int) d_num_blocks.second-1 )
    {
	int i = i_edges.size() - 1;
	int j = j_edges.size() - 1;

	idx_i = my_i_block*d_mesh->getLocalNumCells().first + i;
	idx_j = my_j_block*d_mesh->getLocalNumCells().second + j;
	row_idx = idx_i + idx_j*d_global_edges.first.size();
	d_local_rows.push_back( row_idx );
    }

    // Set the ghost global rows.
    for ( int j = 0; j < (int) j_edges.size(); ++j )
    {
	for ( int i = 0; i < (int) i_edges.size(); ++i )
	{
	    idx_i = my_i_block*d_mesh->getLocalNumCells().first + i;
	    idx_j = my_j_block*d_mesh->getLocalNumCells().second + j;
	    row_idx = idx_i + idx_j*d_global_edges.first.size();
	    d_ghost_local_rows.push_back( row_idx );
	}
    }
    std::cout << comm->getRank() << ": " << d_local_rows << std::endl;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Partitioner::~Partitioner()
{ /* ... */ }

//---------------------------------------------------------------------------//

} 

// end namespace Chimera

//---------------------------------------------------------------------------//
// end Partitioner.cpp
//---------------------------------------------------------------------------//

