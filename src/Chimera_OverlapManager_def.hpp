//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*! 
 * \file Chimera_OverlapManager_def.hpp
 * \author Stuart R. Slattery
 * \brief Overlapping-Domain Decomposition Manager definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_OVERLAPMANAGER_DEF_HPP
#define Chimera_OVERLAPMANAGER_DEF_HPP

#include <algorithm>

#include "Chimera_Assertion.hpp"
#include "Chimera_OperatorTools.hpp"

#include <Teuchos_Array.hpp>

#include <Tpetra_Map.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO>
OverlapManager<Scalar,LO,GO>::OverlapManager( 
    const RCP_TpetraCrsMatrix& iteration_matrix,
    const RCP_TpetraCrsMatrix& probability_matrix,
    const RCP_TpetraVector& lhs,
    const RCP_ParameterList& plist )
    : d_num_overlap( plist->get<GO>("NUM OVERLAP") )
    , d_overlap_iteration_matrix( 
	Tpetra::createCrsMatrix<Scalar,LO,GO>( iteration_matrix->getRowMap() ) )
{
    // Setup for overlap construction.
    Teuchos::ArrayView<const GO> global_rows;
    typename Teuchos::ArrayView<const GO>::const_iterator global_rows_it;
    typename Teuchos::Array<GO>::iterator ghost_global_bound;

    // Get the initial off proc columns.
    Teuchos::Array<GO> ghost_global_rows = 
	OperatorTools::getOffProcColumns( iteration_matrix );

    // Build the overlap first in the iteration matrix by traversing the
    // graph.
    for ( GO i = 0; i < d_num_overlap; ++i )
    {
	// Get rid of the global rows that belong to the original iteration
	// matrix. We don't need to store these, just the overlap.
	global_rows = iteration_matrix->getRowMap()->getNodeElementList();
	for ( global_rows_it = global_rows.begin();
	      global_rows_it != global_rows.end();
	      ++global_rows_it )
	{
	    ghost_global_bound = std::remove( ghost_global_rows.begin(), 
					      ghost_global_rows.end(), 
					      *global_rows_it );
	    ghost_global_rows.resize( std::distance(ghost_global_rows.begin(),
						    ghost_global_bound) );
	}

	// Get the current set of global rows in the overlap probability
	// matrix. 
	global_rows = 
	    d_overlap_iteration_matrix->getRowMap()->getNodeElementList();

	// Append the off proc overlap columns to the global rows.
	for ( global_rows_it = global_rows.begin();
	      global_rows_it != global_rows.end();
	      ++global_rows_it )
	{
	    ghost_global_rows.push_back( *global_rows_it );
	}
	
	// Make a new map of the combined global rows and off proc columns.
	Teuchos::RCP<const Tpetra::Map<LO,GO> > ghost_map = 
	    Tpetra::createNonContigMap<LO,GO>( 
	    ghost_global_rows(), d_overlap_iteration_matrix->getComm() );

	// Export the overlap iteration matrix with the new overlap.
	Tpetra::Export<LO,GO> ghost_exporter( 
	    d_overlap_iteration_matrix->getRowMap(), ghost_map );

	// Update the overlap iteration matrix for the next iteration.
	d_overlap_iteration_matrix = Tpetra::exportAndFillCompleteCrsMatrix<
	    Tpetra::CrsMatrix<Scalar,LO,GO> >(
		iteration_matrix, ghost_exporter );

	// Get the off proc columns.
	ghost_global_rows = 
	    OperatorTools::getOffProcColumns( d_overlap_iteration_matrix );
    }

    // Build the base-to-overlap export.
    Tpetra::Export<LO,GO> base_to_overlap_export( 
	iteration_matrix->getRowMap(), d_overlap_iteration_matrix->getRowMap() );

    // Build the overlap-to-base exporter.
    d_overlap_to_base_export = Teuchos::rcp( 
	new Tpetra::Export<LO,GO>( d_overlap_iteration_matrix->getRowMap(), 
				   iteration_matrix->getRowMap() ) );

    // Apply the overlap to the probability matrix.
    d_overlap_probability_matrix =
	Tpetra::exportAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<Scalar,LO,GO> >(
	    probability_matrix, base_to_overlap_export );

    // Apply the overlap to the LHS.
    d_overlap_lhs = Tpetra::createVector<Scalar,LO,GO>( 
	d_overlap_iteration_matrix->getRowMap() );
    d_overlap_lhs->doExport( *lhs, base_to_overlap_export, Tpetra::INSERT );

    // Move one more step in the graph to get the ghost elements for the
    // iteration matrix.
    global_rows = d_overlap_iteration_matrix->getRowMap()->getNodeElementList();
    ghost_global_rows = 
	OperatorTools::getOffProcColumns( d_overlap_iteration_matrix );

    // Remove the original rows again to only store the overlap data plus one
    // ghosted state in the iteration matrix.
    global_rows = iteration_matrix->getRowMap()->getNodeElementList();
    for ( global_rows_it = global_rows.begin();
	  global_rows_it != global_rows.end();
	  ++global_rows_it )
    {
	ghost_global_bound = std::remove( ghost_global_rows.begin(), 
					  ghost_global_rows.end(), 
					  *global_rows_it );
	ghost_global_rows.resize( std::distance(ghost_global_rows.begin(),
						ghost_global_bound) );
    }

    // Append the off proc overlap columns.
    for ( global_rows_it = global_rows.begin();
	  global_rows_it != global_rows.end();
	  ++global_rows_it )
    {
	ghost_global_rows.push_back( *global_rows_it );
    }

    // Build the new iteration matrix map.
    Teuchos::RCP<const Tpetra::Map<LO,GO> > overlap_ghost_map = 
	Tpetra::createNonContigMap<LO,GO>( 
	    ghost_global_rows(), d_overlap_iteration_matrix->getComm() );

    // Export the iteration matrix to its final form.
    Tpetra::Export<LO,GO> overlap_ghost_exporter( 
	d_overlap_iteration_matrix->getRowMap(), overlap_ghost_map );

    d_overlap_iteration_matrix = Tpetra::exportAndFillCompleteCrsMatrix<
	Tpetra::CrsMatrix<Scalar,LO,GO> >(
	    d_overlap_iteration_matrix, overlap_ghost_exporter );

    // Test postconditions.
    testPostcondition( !d_overlap_iteration_matrix.is_null() );
    testPostcondition( !d_overlap_probability_matrix.is_null() );
    testPostcondition( !d_overlap_lhs.is_null() );
    testPostcondition( !d_overlap_to_base_export.is_null() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO>
OverlapManager<Scalar,LO,GO>::~OverlapManager()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*
 * \brief Determine if a global state is in the overlap owned by this proc. 
 */
template<class Scalar, class LO, class GO>
bool OverlapManager<Scalar,LO,GO>::isOverlapGlobalElement( 
    const GO global_state )
{
    return d_overlap_probability_matrix->getRowMap()->isNodeGlobalElement( 
	global_state );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//

#endif // end Chimera_OVERLAPMANAGER_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_OverlapManager_def.hpp
//---------------------------------------------------------------------------//
