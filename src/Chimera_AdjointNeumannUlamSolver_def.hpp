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
 * \file Chimera_AdjointNeumannUlamSolver_def.hpp
 * \author Stuart R. Slattery
 * \brief Adjoint Neumann-Ulam solver definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_ADJOINTNEUMANNULAMSOLVER_DEF_HPP
#define Chimera_ADJOINTNEUMANNULAMSOLVER_DEF_HPP

#include "Chimera_Assertion.hpp"
#include "Chimera_SamplingTools.hpp"
#include "Chimera_OperatorTools.hpp"
#include "Chimera_RNGTraits.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Export.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO, class RNG>
AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::AdjointNeumannUlamSolver( 
    const RCP_LinearProblem& linear_problem,
    const RCP_LinearOperatorSplit& linear_operator_split,
    const RCP_RNG& rng,
    const RCP_ParameterList& plist )
    : d_relative_weight_cutoff( plist->get<Scalar>("WEIGHT CUTOFF") )
{
    // Build the base class data.
    this->b_linear_problem = linear_problem;
    this->b_linear_operator_split = linear_operator_split;
    this->b_rng = rng;
    this->b_weight_cutoff = d_relative_weight_cutoff;
    this->b_histories_per_stage = plist->get<int>("HISTORIES PER STAGE");

    // Set a unique RNG seed for this process.
    RNGTraits<RNG>::setSeed( 
	*(this->b_rng), 
	this->b_histories_per_stage * 
	this->b_linear_problem->getOperator()->getComm()->getSize() );

    // Build the probability matrix.
    buildProbabilityMatrix();
    testPostcondition( !d_probability_matrix.is_null() );

    // Build the overlap manager.
    d_overlap_manager = Teuchos::rcp( 
	new OverlapManagerType(
	    this->b_linear_operator_split->iterationMatrix(),
	    d_probability_matrix,
	    this->b_linear_problem->getLHS(),
	    plist ) );
    testPostcondition( !d_overlap_manager.is_null() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor
 */
template<class Scalar, class LO, class GO, class RNG>
AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::~AdjointNeumannUlamSolver()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Execute a stage of histories with random walks.
 */
template<class Scalar, class LO, class GO, class RNG>
void AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::walk()
{ 
    // Get the state map and column map.
    Teuchos::RCP<const Tpetra::Map<LO,GO> > state_map =
	this->b_linear_problem->getOperator()->getRowMap();

    Teuchos::RCP<const Tpetra::Map<LO,GO> > col_map =
	this->b_linear_problem->getOperator()->getColMap();

    // Zero out the solution vector.
    this->b_linear_problem->getLHS()->putScalar( 0.0 );

    // Sample the source and populate the history bank.
    HistoryBank<HistoryType> bank = sampleSource();

    // Setup a history buffer.
    HistoryBuffer<HistoryType> buffer;

    // Random walk setup.
    LO local_state = 0;
    LO new_local_state = 0;
    GO global_state = 0;
    GO new_global_state = 0;
    Scalar transition_h, transition_p, transition_weight;
    bool new_state_is_local = false;
    bool state_in_overlap = false;
    bool walk = true;

    // Random walk until all global histories in the stage are terminated.
    while ( walk )
    {
	// If the bank isn't empty, process the top history. We need this here
	// so processes that have emptied their banks can wait until all banks
	// are empty and buffer communication can occur.
	if ( !bank.empty() )
	{
	    // Get the current history state.
	    global_state = bank.top().globalState();
	    local_state = state_map->getLocalElement( global_state );

	    // Check if the current state is in the overlap.
	    state_in_overlap = 
		d_overlap_manager->isOverlapGlobalElement( global_state );

	    // Update LHS tally.
	    if ( state_in_overlap )
	    {
		d_overlap_manager->getOverlapLHS()->sumIntoGlobalValue( 
		    global_state, bank.top().weight() );
	    }
	    else
	    {
		this->b_linear_problem->getLHS()->sumIntoGlobalValue( 
		    global_state, bank.top().weight() );
	    }

	    // Sample the probability matrix to get the new state.
	    new_local_state = 
		sampleProbabilityMatrix( local_state, state_in_overlap );

	    // Invalid state. Terminate the history.
	    if ( new_local_state == Teuchos::OrdinalTraits<LO>::invalid() )
	    {
		transition_weight = 0.0;
	    }
	    // The new state is valid.
	    else
	    {
		new_global_state = col_map->getGlobalElement( new_local_state );

		// Check if the new state is on process.
		new_state_is_local = 
		    state_map->isNodeGlobalElement( new_global_state );

		// Get the iteration matrix transition data.
		if ( new_state_is_local )
		{
		    transition_h = OperatorTools::getMatrixComponentFromLocal(
			this->b_linear_operator_split->iterationMatrix(), 
			new_local_state, local_state );
		}
		else
		{
		    transition_h = OperatorTools::getMatrixComponentFromGlobal(
			d_overlap_manager->getOverlapIterationMatrix(), 
			new_global_state, global_state );
		}

		// Get the transition probability.
		if ( state_in_overlap )
		{
		    transition_p = OperatorTools::getMatrixComponentFromLocal(
			d_overlap_manager->getOverlapProbabilityMatrix(), 
			local_state, new_local_state );
		}
		else
		{
		    transition_p = OperatorTools::getMatrixComponentFromLocal(
			d_probability_matrix, local_state, new_local_state );
		}

		// Compute the state transition weight.
		if ( transition_p == 0.0 )
		{
		    transition_weight = 0.0;
		}
		else
		{
		    transition_weight = std::abs( transition_h / transition_p );
		}
	    }

	    // Update the history for the transition.
	    bank.top().setGlobalState( new_global_state );
	    bank.top().multiplyWeight( transition_weight );

	    // If the history is below the weight cutoff it is terminated.
	    if ( d_relative_weight_cutoff > bank.top().weightAbs() )
	    {
		bank.pop();
	    }
	    // Else if the history has left the local domain, buffer it.
	    else if ( !new_state_is_local && !state_in_overlap )
	    {
		buffer.pushBack( bank.pop() );
	    }
	}

	// Check if the banks are empty.
	if ( allBanksEmpty( bank ) )
	{
	    // If all the buffers are empty, the stage is complete.
	    if ( allBuffersEmpty( buffer ) )
	    {
		walk = false;
	    }
	    // Else the buffers are not empty, communicate them.
	    else
	    {
		bank = buffer.communicate( state_map );
	    }
	}
    }

    // Export the overlap LHS to the base decomposition LHS. Sum the tallies.
    d_overlap_manager->exportOverlapLHS();

    // Scale the solution by the number of histories in this stage.
    Scalar solution_scaling = 
	1.0 / Teuchos::as<Scalar>(this->b_histories_per_stage);
    this->b_linear_problem->getLHS()->scale( solution_scaling );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the adjoint probability matrix.
 */
template<class Scalar, class LO, class GO, class RNG>
void AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::buildProbabilityMatrix()
{
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> > iteration_matrix =
	this->b_linear_operator_split->iterationMatrix();
    Teuchos::RCP<const Tpetra::Map<LO,GO> > row_map = 
	iteration_matrix->getRowMap();
    Teuchos::RCP<const Tpetra::Map<LO,GO> > col_map = 
	iteration_matrix->getColMap();
    d_probability_matrix = Tpetra::createCrsMatrix<Scalar>( row_map );

    Scalar row_sum = 0.0;

    Teuchos::ArrayView<const LO> row_indices;
    typename Teuchos::ArrayView<const LO>::const_iterator row_indices_it;
    Teuchos::ArrayView<const Scalar> row_values;
    typename Teuchos::ArrayView<const Scalar>::const_iterator row_values_it;

    Teuchos::Array<GO> probability_col(1);
    Teuchos::Array<Scalar> probability_value(1);

    LO local_row;
    GO probability_row;
    Teuchos::ArrayView<const GO> global_rows = row_map->getNodeElementList();
    typename Teuchos::ArrayView<const GO>::const_iterator row_it;
    for ( row_it = global_rows.begin(); row_it != global_rows.end(); ++row_it )
    {
	local_row = row_map->getLocalElement( *row_it );

	iteration_matrix->getLocalRowView( 
	    local_row, row_indices, row_values );

	row_sum = 0.0;
	for ( row_values_it = row_values.begin();
	      row_values_it != row_values.end();
	      ++row_values_it )
	{
	    row_sum += std::abs( *row_values_it );
	}

	probability_col[0] = *row_it;

	for ( row_indices_it = row_indices.begin(),
	       row_values_it = row_values.begin();
	      row_indices_it != row_indices.end();
	      ++row_indices_it, ++row_values_it )
	{
	    if ( row_sum > 0.0 )
	    {
		probability_value[0] = std::abs( *row_values_it ) / row_sum;
	    }
	    else
	    {
		probability_value[0] = 0.0;
	    }

	    probability_row = col_map->getGlobalElement( *row_indices_it );

	    d_probability_matrix->insertGlobalValues( probability_row,
						      probability_col(),
						      probability_value() );
	}
    }

    d_probability_matrix->fillComplete();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample the source to build a starting history bank and set the
 * random number generator seeds.
 */
template<class Scalar, class LO, class GO, class RNG>
HistoryBank<History<Scalar,GO> > 
AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::sampleSource()
{
    // Get the RHS.
    Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> > source = 
	this->b_linear_problem->getRHS();
    Teuchos::RCP<const Tpetra::Map<LO,GO> > source_map = source->getMap();
    Teuchos::ArrayView<const GO> global_states = 
	source_map->getNodeElementList();
    Teuchos::ArrayRCP<const Scalar> local_source_view =
	source->get1dView();

    // Stratify sample the source.
    Teuchos::ArrayRCP<GO> starting_states = 
	SamplingTools::stratifySampleGlobalPDF( 
	    this->b_histories_per_stage, source );
    testInvariant( local_source_view.size() == starting_states.size() );

    // Get the starting source weight.
    Scalar source_weight = source->norm1();

    // Set the relative weight cutoff.
    d_relative_weight_cutoff *= source_weight;

    // Build the starting source bank.
    HistoryBank<HistoryType> source_bank;
    Scalar history_weight = 0.0;
    typename Teuchos::ArrayRCP<const Scalar>::const_iterator local_source_it;
    typename Teuchos::ArrayRCP<GO>::const_iterator starting_states_it;
    typename Teuchos::ArrayView<const GO>::const_iterator global_states_it;
    for ( local_source_it = local_source_view.begin(),
       starting_states_it = starting_states.begin(),
	 global_states_it = global_states.begin();
	  local_source_it != local_source_view.end();
	  ++local_source_it, ++starting_states_it, ++global_states_it )
    {
	history_weight = source_weight * (*local_source_it) 
			 / std::abs(*local_source_it);

	for ( GO i = 0; i < *starting_states_it; ++i )
	{
	    source_bank.push( 
		HistoryType( history_weight, *global_states_it ) );
	}
    }

    return source_bank;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample a row of the probability matrix to get a new local state.
 */
template<class Scalar, class LO, class GO, class RNG>
LO AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::sampleProbabilityMatrix(
    const LO local_state, const bool state_in_overlap )
{
    Teuchos::ArrayView<const LO> local_indices;
    Teuchos::ArrayView<const Scalar> local_values;

    if ( state_in_overlap )
    {
	d_overlap_manager->getOverlapProbabilityMatrix()->getLocalRowView(
	    local_state, local_indices, local_values );
    }
    else
    {
	d_probability_matrix->getLocalRowView( 
	    local_state, local_indices, local_values );
    }

    return SamplingTools::sampleLocalDiscretePDF(
	local_values, local_indices, this->b_rng );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check for completion of all local random walks in the bank.
 */
template<class Scalar, class LO, class GO, class RNG>
bool AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::allBanksEmpty(
    const HistoryBank<HistoryType>& bank )
{
    typedef typename HistoryBank<HistoryType>::size_type size_type;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	this->b_linear_problem->getOperator()->getComm();

    size_type global_num_histories = 0;

    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, bank.size(), 
			Teuchos::Ptr<size_type>(&global_num_histories) );

    return 0 == global_num_histories;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check for empty buffers on all processes.
 */
template<class Scalar, class LO, class GO, class RNG>
bool AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::allBuffersEmpty(
    const HistoryBuffer<HistoryType>& buffer )
{
    typedef typename HistoryBuffer<HistoryType>::size_type size_type;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	this->b_linear_problem->getOperator()->getComm();

    size_type global_num_histories = 0;

    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, buffer.size(), 
			Teuchos::Ptr<size_type>(&global_num_histories) );

    return 0 == global_num_histories;
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_ADJOINTNEUMANNULAMSOLVER_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_AdjointNeumannUlamSolver_def.hpp
//---------------------------------------------------------------------------//
