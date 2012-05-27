//---------------------------------------------------------------------------//
/*!
 * \file Mesh.cpp
 * \author Stuart R. Slattery
 * \brief 2D Cartesian mesh definition.
 */
//---------------------------------------------------------------------------//

#include <cassert>

#include "Mesh.hpp"

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Mesh::Mesh( long_type global_Ni, long_type global_Nj,
	    double global_i_min, double global_j_min, 
	    double global_i_max, double global_j_max,
	    size_type I_block, size_type J_block,
	    std::vector<double> i_edges, std::vector<double> j_edges )
    : d_global_N( global_Ni, global_Nj )
    , d_lower( global_i_min, global_j_min )
    , d_upper( global_i_max, global_j_max )
    , d_block( I_block, J_block )
    , d_local_N( i_edges.size(), j_edges.size() )
    , d_edges( i_edges, j_edges )
{
    computeCellWidths();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
Mesh::~Mesh()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Compute cell widths.
 */
void Mesh::computeCellWidths()
{
    std::vector<double>::const_iterator edge_it;

    for ( edge_it = d_edges.first.begin(); 
	  edge_it != (d_edges.first.end()-1); 
	  ++edge_it )
    {
	d_widths.first.push_back( *(edge_it+1) - *(edge_it) );
    }
    
    assert( d_widths.first.size() == d_edges.first.size() - 1 );

    for ( edge_it = d_edges.second.begin(); 
	  edge_it != (d_edges.second.end()-1); 
	  ++edge_it )
    {
	d_widths.second.push_back( *(edge_it+1) - *(edge_it) );
    }

    assert( d_widths.second.size() == d_edges.second.size() - 1 );
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Mesh.cpp
//---------------------------------------------------------------------------//

