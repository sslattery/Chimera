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
 * \file Chimera_AdjointNeumannUlamSolver.hpp
 * \author Stuart R. Slattery
 * \brief Adjoint Neumann-Ulam solver definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_ADJOINTNEUMANNULAMSOLVER_DEF_HPP
#define Chimera_ADJOINTNEUMANNULAMSOLVER_DEF_HPP

#include <numeric>

#include "Chimera_Assertion.hpp"
#include "Chimera_SamplingTools.hpp"
#include "Chimera_OperatorTools.hpp"

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO, class RNG>
AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::AdjointNeumannUlamSolver( 
    const RCP_LinearProblem& linear_problem,
    const RCP_LinearOperatorSplit& operator_split,
    const RCP_RNG& rng,
    const RCP_Parameterlist& plist )
    : d_relative_weight_cutoff( plist->get<Scalar>("WEIGHT CUTOFF") )
{
    this->b_linear_problem = linear_problem;
    this->b_operator_split = operator_split;
    this->b_rng = rng;
    this->b_weight_cutoff = d_relative_weight_cutoff;
    this->b_histories_per_stage = plist->get<int>("HISTORIES PER STAGE");

    buildProbabilityMatrix();
    testPostcondition( !d_probability_matrix.is_null() );
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
    // Get the state map.
    Teuchos::RCP<const Tpetra::Map<LO,GO> > state_map =
	this->b_linear_problem->getOperator->getRowMap();

    // Zero out the local solution vector and get a view to write into.
    this->b_linear_problem->getLHS()->putScalar( 0.0 );
    Teuchos::ArrayRCP<Scalar> lhs_view =     
	this->b_linear_problem->getLHS()->get1dViewNonConst();
   
    // Sample the source and populate the history bank.
    HistoryBank<HistoryType> bank = sampleSource();

    // Setup a history buffer.
    HistoryBuffer<HistoryType> buffer;

    // Random walk until all global histories are terminated.
    LO local_state, new_local_state;
    GO global_state, new_global_state;
    Scalar h_transition, p_transition, transition_weight;
    Teuchos::ArrayView<LO> local_indices;
    Teuchos::ArrayView<Scalar> local_values;
    bool walk = true;
    while ( walk )
    {
	// Get the current state.
	global_state = bank.top().globalState();
	local_state = state_map->getLocalElement( global_state );
	
	// Update LHS.
	lhs_view[ local_state ] += bank.top().weight();

	// Sample the probability matrix to get the new state.
	this->b_linear_problem->getOperator()->localRowView( 
	    local_state, local_indices, local_values );
	new_local_state = SamplingTools::sampleLocalDiscretePDF(
	    local_values, local_indices, d_rng );

	// Compute state transition weight.
	new_global_state = state_map->getGlobalElement( new_local_state );
	h_transition = OperatorTools::getMatrixComponentFromLocal(
	    this->b_linear_problem->getOperator(), 
	    new_local_state, local_state );
	p_transition = OperatorTools::getMatrixComponentsFromLocal(
	    d_probability_matrix, local_state, new_local_state );
	transition_weight = h_transition / p_transition;

	// Update the history for the transition.
	bank.top().setGlobalState( new_global_state );
	bank.top().multiplyWeight( transition_weight );

	// If the history is below the weight cutoff it is terminated.
	if ( d_relative_weight_cutoff > bank.top().weightAbs() )
	{
	    bank.pop();
	}

	// If the history has left the local domain, buffer it.
	else if ( !state_map->isNodeGlobalElement( bank.top().globalState() ) )
	{
	    buffer.push_back( bank.pop() );
	}

	// Check if the banks are empty.
	if ( allRandomWalksComplete( bank ) )
	{
	    // If they are empty and the all the buffers are empty, the walk
	    // is complete.
	    if ( allBuffersEmpty( buffer ) )
	    {
		walk = false;
	    }

	    // If the buffers are not empty, communicate them.
	    else
	    {
		bank = buffer.communicate( state_map );
	    }
	}
    }

    // Scale the solution by the number of histories in this stage.
    Scalar solution_scaling = 1.0 / this->b_histories_per_stage;
    this->b_linear_problem->getLHS()->scale( solution_scaling );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the probability matrix.
 */
template<class Scalar, class LO, class GO, class RNG>
void AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::buildProbabilityMatrix()
{
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> iteration_matrix =
	this->b_linear_operator_split->iterationMatrix();
    Teuchos::RCP<const Tpetra::Map<LO,GO> > row_map = 
	iteration_matrix->getRowMap();
    d_probability_matrix = Tpetra::createCrsMatrix<Scalar>( row_map );

    Scalar row_sum = 0.0;

    Teuchos::ArrayView<const GO> row_indices;
    typename Teuchos::ArrayView<const GO>::const_iterator row_indices_it;
    Teuchos::ArrayView<const Scalar> row_values;
    typename Teuchos::ArrayView<const Scalar>::const_iterator row_values_it;

    Teuchos::Array<GO> probability_col(1);
    Teuchos::Array<Scalar> probability_value(1);

    Teuchos::ArrayView<const GO> global_rows = row_map->getNodeElementList();
    typename Teuchos::ArrayView<const GO>::const_iterator row_it;
    for ( row_it = global_rows.begin(); row_it != global_rows.end(); ++row_it )
    {
	iteration_matrix->getGlobalRowView( 
	    *row_it, row_indices, row_values );

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

	    d_probability_matrix->insertGlobalValues( *row_indices_it,
						      probability_col(),
						      probability_value() );
	}
    }

    d_probability_matrix->fillComplete();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Sample the source to build a starting history bank.
 */
template<class Scalar, class LO, class GO, class RNG>
HistoryBank<History<Scalar,GO> > 
AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::sampleSource()
{
    Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> > source = 
	this->b_linear_problem->getRHS();

    Teuchos::ArrayRCP<Scalar> local_source_view = source->get1dView();

    Teuchos::RCP<const Tpetra::Map<LO,GO> > source_map = source->getMap();
    Teuchos::ArrayView<const GO> global_states = 
	source_map->getNodeElementList();
    testInvariant( local_source_view.size() == global_states.size() );

    Teuchos::ArrayRCP<GO> starting_states = 
	SamplingTools::stratifySampleGlobalPDF( this->b_histories_per_stage,
						source );
    testInvariant( local_source_view.size() == starting_states.size() );

    HistoryBank<HistoryType> source_bank;
    Scalar source_weight = source->norm1();
    d_relative_weight_cutoff *= source_weight;

    Scalar history_weight = 0.0;
    typename Teuchos::ArrayRCP<Scalar>::const_iterator local_source_it;
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
 * \brief Check for completion of all random walks.
 */
template<class Scalar, class LO, class GO, class RNG>
bool AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::allRandomWalksComplete(
    const HistoryBank<HistoryType>& bank )
{
    typedef typename HistoryBank<HistoryType>::size_type size_type;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	this->b_linear_problem->getOperator()->getComm();

    size_type global_num_histories = 0;

    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, bank.size(), 
			Teuchos::Ptr<size_type>(*global_num_histories) );

    return 0 == global_num_histories;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Check for empty buffers on all processes.
 */
template<class Scalar, class LO, class GO, class RNG>
bool AdjointNeumannUlamSolver<Scalar,LO,GO,RNG>::allBuffersEmpty
    const HistoryBuffer<HistoryType>& buffer )
{
    typedef typename HistoryBuffer<HistoryType>::size_type size_type;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	this->b_linear_problem->getOperator()->getComm();

    size_type global_num_histories = 0;

    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, buffer.size(), 
			Teuchos::Ptr<size_type>(*global_num_histories) );

    return 0 == global_num_histories;
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_ADJOINTNEUMANNULAMSOLVER_HPP

//---------------------------------------------------------------------------//
// end Chimera_AdjointNeumannUlamSolver.hpp
//---------------------------------------------------------------------------//
