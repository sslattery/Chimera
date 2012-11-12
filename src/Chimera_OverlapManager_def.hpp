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
    : d_lhs( lhs )
    , d_num_overlap( plist->get<GO>("NUM OVERLAP") )
{
    testPrecondition( d_num_overlap >= 0 );

    // Create the initial empty overlap iteration matrix.
    Teuchos::RCP<const Tpetra::Map<int> > empty_map = 
	Tpetra::createUniformContigMap<int,int>( 
	    0, iteration_matrix->getComm() );
    d_overlap_iteration_matrix = 
	Tpetra::createCrsMatrix<Scalar,LO,GO>( empty_map );
    d_overlap_iteration_matrix->fillComplete();

    // Build the overlap if number of overlapping states is greater than 0.
    if ( d_num_overlap > 0 )
    {
	buildOverlap( iteration_matrix, probability_matrix );
    }

    // Build the ghost rows for the iteration matrix.
    buildIterationGhost( iteration_matrix );
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
 * \brief Export the overlap LHS to the original LHS by summing the tallies.
 */
template<class Scalar, class LO, class GO>
void OverlapManager<Scalar,LO,GO>::exportOverlapLHS()
{
    if ( d_num_overlap > 0 )
    {
	d_lhs->doExport( 
	    *d_overlap_lhs, *d_overlap_to_base_export, Tpetra::ADD );
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Determine if a global state is in the overlap owned by this proc. 
 */
template<class Scalar, class LO, class GO>
bool OverlapManager<Scalar,LO,GO>::isOverlapGlobalElement( 
    const GO global_state )
{
    if ( d_num_overlap > 0 )
    {
	return d_overlap_probability_matrix->getRowMap()->isNodeGlobalElement(
	    global_state );
    }
    else
    {
	return false;
    }
}

//---------------------------------------------------------------------------//
/*
 * \brief Build the overlap.
 */
template<class Scalar, class LO, class GO>
void OverlapManager<Scalar,LO,GO>::buildOverlap( 
    const RCP_TpetraCrsMatrix& iteration_matrix,
    const RCP_TpetraCrsMatrix& probability_matrix )
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

	// Append the on proc overlap columns to the off proc columns.
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
	    iteration_matrix->getRowMap(), ghost_map );

	// Update the overlap iteration matrix for the next iteration.
	d_overlap_iteration_matrix = Tpetra::exportAndFillCompleteCrsMatrix<
	    Tpetra::CrsMatrix<Scalar,LO,GO> >(
		iteration_matrix, ghost_exporter );

	// Get the next rows in the graph.
	ghost_global_rows = 
	    OperatorTools::getOffProcColumns( d_overlap_iteration_matrix );
    }

    // Apply the overlap to the probability matrix.
    Tpetra::Export<LO,GO> probability_export( 
	probability_matrix->getRowMap(), d_overlap_iteration_matrix->getRowMap() );
    d_overlap_probability_matrix =
	Tpetra::exportAndFillCompleteCrsMatrix<Tpetra::CrsMatrix<Scalar,LO,GO> >(
	    probability_matrix, probability_export );

    // Apply the overlap to the LHS.
    d_overlap_lhs = Tpetra::createVector<Scalar,LO,GO>( 
	d_overlap_iteration_matrix->getRowMap() );
    Tpetra::Export<LO,GO> lhs_export( d_lhs->getMap(), d_overlap_lhs->getMap() );
    d_overlap_lhs->doExport( *d_lhs, lhs_export, Tpetra::INSERT );

    // Build the overlap-to-base exporter for the LHS.
    d_overlap_to_base_export = Teuchos::rcp( 
	new Tpetra::Export<LO,GO>( d_overlap_lhs->getMap(), d_lhs->getMap() ) );

    // Test postconditions.
    testPostcondition( !d_overlap_iteration_matrix.is_null() );
    testPostcondition( !d_overlap_probability_matrix.is_null() );
    testPostcondition( !d_overlap_lhs.is_null() );
    testPostcondition( !d_overlap_to_base_export.is_null() );
}

//---------------------------------------------------------------------------//
/*
 * \brief Build the iteration ghost rows.
 */
template<class Scalar, class LO, class GO>
void OverlapManager<Scalar,LO,GO>::buildIterationGhost( 
    const RCP_TpetraCrsMatrix& iteration_matrix )
{
    // Setup.
    Teuchos::ArrayView<const GO> global_rows;
    typename Teuchos::ArrayView<const GO>::const_iterator global_rows_it;
    typename Teuchos::Array<GO>::iterator ghost_global_bound;
    Teuchos::Array<GO> ghost_global_rows;

    // Move one more step in the graph to get the ghost elements for the
    // iteration matrix.
    if ( d_num_overlap == 0 )
    {
	ghost_global_rows = 
	    OperatorTools::getOffProcColumns( iteration_matrix );
    }
    else
    {
	// Get the next rows in the graph.
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

	// Append the on proc overlap columns.
	global_rows = 
	    d_overlap_iteration_matrix->getRowMap()->getNodeElementList();
	for ( global_rows_it = global_rows.begin();
	      global_rows_it != global_rows.end();
	      ++global_rows_it )
	{
	    ghost_global_rows.push_back( *global_rows_it );
	}
    }

    // Build the new iteration matrix map.
    Teuchos::RCP<const Tpetra::Map<LO,GO> > overlap_ghost_map = 
	Tpetra::createNonContigMap<LO,GO>( 
	    ghost_global_rows(), d_overlap_iteration_matrix->getComm() );

    // Export the iteration matrix to its final form.
    Tpetra::Export<LO,GO> overlap_ghost_exporter( 
	iteration_matrix->getRowMap(), overlap_ghost_map );

    d_overlap_iteration_matrix = Tpetra::exportAndFillCompleteCrsMatrix<
	Tpetra::CrsMatrix<Scalar,LO,GO> >(
	    iteration_matrix, overlap_ghost_exporter );

    testPostcondition( !d_overlap_iteration_matrix.is_null() );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//

#endif // end Chimera_OVERLAPMANAGER_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_OverlapManager_def.hpp
//---------------------------------------------------------------------------//
