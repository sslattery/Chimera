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

#include "Partitioner.hpp"

#include <Exception.hpp>

#include <Teuchos_ENull.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Partitioner::Partitioner( RCP_Comm comm, RCP_ParameterList plist )
    : d_num_blocks( SizePair( plist->get<int>( "I_BLOCKS" ), 
			      plist->get<int>( "J_BLOCKS" ) ) )
{
    // Comm parameters.
    unsigned int myRank = comm->getRank();
    unsigned int mySize = comm->getSize();

    // Check that the block specification and communicator are consistent.
    testPrecondition( d_num_blocks.first * d_num_blocks.second == mySize,
		      "I_BLOCKS*J_BLOCKS != number of processors" );

    // Block indices.
    int myJBlock = std::floor( myRank / d_num_blocks.first );
    int myIBlock = myRank - d_num_blocks.first*myJBlock;

    std::vector<double> i_edges, j_edges;
    double global_i_min, global_i_max, global_j_min, global_j_max;
    int global_num_i, global_num_j;

    // Uniform grid case.
    if ( plist->get<std::string>("GRID_TYPE") == "UNIFORM" )
    {
	global_i_min = plist->get<double>( "X_MIN" );
	global_i_max = plist->get<double>( "X_MAX" );
	global_j_min = plist->get<double>( "Y_MIN" );
	global_j_max = plist->get<double>( "Y_MAX" );
	global_num_i = plist->get<int>( "X_NUM_CELLS" );
	global_num_j = plist->get<int>( "Y_NUM_CELLS" );

	double width_i = (global_i_max - global_i_min) / global_num_i;
	double width_j = (global_j_max - global_j_min) / global_num_j;

	int i_edges_size = std::floor( global_num_i / d_num_blocks.first ) + 1;
	int j_edges_size = std::floor( global_num_j / d_num_blocks.second ) + 1;

	int i_remainder = global_num_i % d_num_blocks.first;
	int j_remainder = global_num_j % d_num_blocks.second;

	double i_edge_val = 0.0;
	for ( int i = 0; i < i_edges_size; ++i )
	{
	    i_edges.push_back( myRank*(i_edges_size-1) + i_edge_val );
	    i_edge_val += width_i;
	}
	if ( myRank == mySize - 1 )
	{
	    for ( int i = 0; i < i_remainder; ++i )
	    {
		i_edges.push_back( myRank*(i_edges_size-1) + i_edge_val );
		i_edge_val += width_i;
	    }
	}
	comm->barrier();

	double j_edge_val = 0.0;
	for ( int j = 0; j < j_edges_size; ++j )
	{
	    j_edges.push_back( myRank*(j_edges_size-1) + j_edge_val );
	    j_edge_val += width_j;
	}
	if ( myRank == mySize - 1 )
	{
	    for ( int j = 0; j < j_remainder; ++j )
	    {
		j_edges.push_back( myRank*(j_edges_size-1) + j_edge_val );
		j_edge_val += width_j;
	    }
	}
	comm->barrier();
    }
    
    // Nonuniform grid case.
    else if ( plist->get<std::string>("GRID_TYPE") == "NONUNIFORM" )
    {
	std::vector<double> global_i_edges = 
	    plist->get< std::vector<double> >("I_EDGES");
	std::vector<double> global_j_edges 
	    = plist->get< std::vector<double> >("J_EDGES");

	global_i_min = global_i_edges.front();
	global_i_max = global_i_edges.back();
	global_j_min = global_j_edges.front();
	global_j_max = global_j_edges.back();
	global_num_i = global_i_edges.size() - 1;
	global_num_j = global_j_edges.size() - 1;

	int i_edges_size = std::floor( global_num_i / d_num_blocks.first ) + 1;
	int j_edges_size = std::floor( global_num_j / d_num_blocks.second ) + 1;

	int i_remainder = global_num_i % d_num_blocks.first;
	int j_remainder = global_num_j % d_num_blocks.second;

	int i_edge_idx;
	for ( int i = 0; i < i_edges_size; ++i )
	{
	    i_edges.push_back( i_edges[ myRank*(i_edges_size-1) + i_edge_idx ] );
	    ++i_edge_idx;
	}
	if ( myRank == mySize - 1 )
	{
	    for ( int i = 0; i < i_remainder; ++i )
	    {
		i_edges.push_back( i_edges[ myRank*(i_edges_size-1) + i_edge_idx ] );
		++i_edge_idx;
	    }
	}
	comm->barrier();

	int j_edge_idx;
	for ( int j = 0; j < j_edges_size; ++j )
	{
	    j_edges.push_back( j_edges[ myRank*(j_edges_size-1) + j_edge_idx ] );
	    ++j_edge_idx;
	}
	if ( myRank == mySize - 1 )
	{
	    for ( int j = 0; j < j_remainder; ++j )
	    {
		j_edges.push_back( j_edges[ myRank*(j_edges_size-1) + j_edge_idx ] );
		++j_edge_idx;
	    }
	}
	comm->barrier();
    }

    // Unsupported cases.
    else
    {
	testPrecondition( plist->get<std::string>("GRID_TYPE") == "NONUNIFORM" ||
			  plist->get<std::string>("GRID_TYPE") == "UNIFORM",
			  "Grid type not supported" );
    }

    // Create th mesh.
    d_mesh = Teuchos::rcp( new Mesh( global_num_i, global_num_j,
				     global_i_min, global_j_min,
				     global_i_max, global_j_max,
				     myIBlock, myJBlock,
				     i_edges, j_edges ) );

    testPostcondition( d_mesh != Teuchos::null,
		       "Error partitioning mesh." );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Partitioner::~Partitioner()
{ /* ... */ }

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Partitioner.cpp
//---------------------------------------------------------------------------//

