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
// \file Chimera_AdjointMC.cpp
// \author Stuart Slattery
// \brief Adjoint Monte Carlo solver definition.
//---------------------------------------------------------------------------//

#include <cmath>
#include <cassert>
#include <vector>
#include <iterator>
#include <algorithm>

#include "Chimera_AdjointMCDev.hpp"

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
AdjointMC::AdjointMC( Teuchos::RCP<Epetra_LinearProblem> &linear_problem,
		      Teuchos::RCP<Teuchos::ParameterList> &plist )
    : d_linear_problem( linear_problem )
    , d_plist( plist )
    , d_rng( RNGTraits<boost::mt11213b>::create() )
    , d_H( buildH() )
    , d_Q( buildQ() )
    , d_C( buildC() )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
AdjointMC::~AdjointMC()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*! 
 * \brief Solve.
 */
void AdjointMC::walk()
{
    // Get the solver parameters.
    int num_histories = d_plist->get<int>("NUM HISTORIES");
    double weight_cutoff = d_plist->get<double>("WEIGHT CUTOFF");

    // Get the LHS and source.
    Epetra_Vector *x = 
	dynamic_cast<Epetra_Vector*>( d_linear_problem->GetLHS() );
    const Epetra_Vector *b = 
	dynamic_cast<Epetra_Vector*>( d_linear_problem->GetRHS() );
    int N = x->GlobalLength();
    int n_H = d_H->GlobalMaxNumEntries();
    int n_Q = d_Q.GlobalMaxNumEntries();
    int n_C = d_C.GlobalMaxNumEntries();

    // Setup.
    int state;
    int new_state;
    int init_state;
    int new_index;
    double weight;
    double zeta;
    double relative_cutoff;
    bool walk;
    std::vector<double> b_cdf( N );
    std::vector<double> H_values( n_H );
    std::vector<int> H_indices( n_H );
    int H_size;
    std::vector<double> Q_values( n_Q );
    std::vector<int> Q_indices( n_Q );
    int Q_size;
    std::vector<double> C_values( n_C );
    std::vector<int> C_indices( n_C );
    int C_size;
    std::vector<int>::iterator Q_it;
    std::vector<int>::iterator H_it;

    // Build source cdf.
    b_cdf[0] = std::abs((*b)[0]);
    double b_norm = b_cdf[0];
    for ( int i = 1; i < N; ++i )
    {
	b_norm += std::abs((*b)[i]);
	b_cdf[i] = std::abs((*b)[i]) + b_cdf[i-1];	
    }
    for ( int i = 0; i < N; ++i )
    {
	b_cdf[i] /= b_norm;
    }

    // Do random walks for specified number of histories.
    double transitions_per_history = 0.0;
    int max_transitions_in_history = 0;
    int transitions = 0;
    for ( int n = 0; n < num_histories; ++n )
    {
	// Sample the source to get the initial state.
	zeta = (double) RNGTraits<boost::mt11213b>::generate(*d_rng) / 
	       RNGTraits<boost::mt11213b>::max(*d_rng);

	// Line source.
	if ( d_plist->get<bool>("LINE SOURCE") )
	{
	    init_state = d_plist->get<int>("SOURCE STATE");
	}
	// Right-hand side PDF source.
	else
	{
	    init_state = std::distance( 
		b_cdf.begin(),
		std::lower_bound( b_cdf.begin(), b_cdf.end(), zeta ) );
	}

	// Random walk.
	weight = b_norm / std::abs((*b)[init_state]);
	relative_cutoff = weight_cutoff*weight;
	// weight = std::abs((*b)[init_state]);
	// relative_cutoff = weight_cutoff*std::abs((*b)[init_state]);
	state = init_state;
	walk = true;
	transitions = 0;

	if ( d_plist->get<bool>("HISTORY DIAGNOSTICS") )
	{
	    std::cout << "NEW PART" << std::endl;
	}

	while ( walk )
	{
	    // Update LHS.
	    (*x)[state] += (*b)[init_state] * weight / num_histories;

	    // Sample the CDF to get the next state.
	    d_C.ExtractGlobalRowCopy( state, 
	    			      n_C, 
	    			      C_size, 
	    			      &C_values[0], 
	    			      &C_indices[0] );

	    zeta = (double) RNGTraits<boost::mt11213b>::generate(*d_rng) / 
	    	   RNGTraits<boost::mt11213b>::max(*d_rng);

	    new_index = std::distance( 
	    	C_values.begin(),
	    	std::lower_bound( C_values.begin(), 
	    			  C_values.begin()+C_size,
	    			  zeta ) );
	    new_state = C_indices[ new_index ];
	    
	    // Get the state components of Q and H.
	    d_Q.ExtractGlobalRowCopy( state, 
				      n_Q, 
				      Q_size, 
				      &Q_values[0], 
				      &Q_indices[0] );

	    d_H->ExtractGlobalRowCopy( new_state, 
				       n_H, 
				       H_size, 
				       &H_values[0], 
				       &H_indices[0] );

	    Q_it = std::find( Q_indices.begin(),
			      Q_indices.begin()+Q_size,
			      new_state );

	    H_it = std::find( H_indices.begin(),
			      H_indices.begin()+H_size,
			      state );

	    // Compute new weight.
	    if ( Q_values[std::distance(Q_indices.begin(),Q_it)] == 0 ||
		 Q_it == Q_indices.begin()+Q_size ||
		 H_it == H_indices.begin()+H_size )
	    {
		weight = 0.0;
	    }
	    else
	    {
		weight *= H_values[std::distance(H_indices.begin(),H_it)] / 
			  Q_values[std::distance(Q_indices.begin(),Q_it)];
	    }

	    // Check the new weight against the cutoff.
	    if ( weight < relative_cutoff )
	    {
		walk = false;
	    }

	    if ( d_plist->get<bool>("HISTORY DIAGNOSTICS") )
	    {
		std::cout << state << " " << new_state << " " << init_state << " " 
			  << weight << " " << relative_cutoff << " " 
			  << transitions << " " 
			  << H_values[std::distance(H_indices.begin(),H_it)] << " "  
			  << Q_values[std::distance(Q_indices.begin(),Q_it)] 
			  << std::endl;
	    }

	    // Update the state.
	    state = new_state;
	    ++transitions;
	}
	
	// Update diagnostics
	max_transitions_in_history = 
	    std::max( transitions, max_transitions_in_history );
	transitions_per_history += transitions;
    }

    transitions_per_history /= num_histories;

    if ( d_plist->get<bool>("ITERATION DIAGNOSTICS") )
    {
	std::cout << "----------------------------------------------" << std::endl;
	std::cout << "Average transitions per history: " 
		  << transitions_per_history << std::endl;
	std::cout << "Max transitions in a history: "
		  << max_transitions_in_history << std::endl;
	std::cout << "----------------------------------------------" << std::endl;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the iteration matrix.
 */
Teuchos::RCP<Epetra_CrsMatrix> AdjointMC::buildH()
{
    const Epetra_CrsMatrix *A = 
	dynamic_cast<Epetra_CrsMatrix*>( d_linear_problem->GetMatrix() );
    Teuchos::RCP<Epetra_CrsMatrix> H = Teuchos::rcp(
	new Epetra_CrsMatrix( Copy, A->RowMap(), A->GlobalMaxNumEntries() ) );
    int N = A->NumGlobalRows();
    int n_A = A->GlobalMaxNumEntries();
    std::vector<double> A_values( n_A );
    std::vector<int> A_indices( n_A );
    int A_size = 0;
    double local_H;
    bool found_diag = false;
    for ( int i = 0; i < N; ++i )
    {
	A->ExtractGlobalRowCopy( i,
				 n_A, 
				 A_size, 
				 &A_values[0], 
				 &A_indices[0] );

	for ( int j = 0; j < A_size; ++j )
	{
	    if ( i == A_indices[j] )
	    {
		local_H = 1.0 - A_values[j];
		H->InsertGlobalValues( i, 1, &local_H, &A_indices[j] );
		found_diag = true;
	    }
	    else
	    {
		local_H = -A_values[j];
		H->InsertGlobalValues( i, 1, &local_H, &A_indices[j] );
	    }

	    if ( !found_diag )
	    {
		local_H = 1.0;
		H->InsertGlobalValues( i, 1, &local_H, &i );
	    }
	}
    }

    H->FillComplete();
    H->OptimizeStorage();
    return H;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the adjoint probability matrix.
 */
Epetra_CrsMatrix AdjointMC::buildQ()
{
    Epetra_CrsMatrix Q( Copy, d_H->RowMap(), d_H->GlobalMaxNumEntries() );
    int N = d_H->NumGlobalRows();
    int n_H = d_H->GlobalMaxNumEntries();
    std::vector<double> H_values( n_H );
    std::vector<int> H_indices( n_H );
    std::vector<int>::iterator H_it;
    int H_size = 0;
    double local_Q = 0.0;
    double row_sum = 0.0;
    for ( int i = 0; i < N; ++i )
    {
	d_H->ExtractGlobalRowCopy( i,
				   n_H, 
				   H_size, 
				   &H_values[0], 
				   &H_indices[0] );

	row_sum = 0.0;
	for ( int j = 0; j < H_size; ++j )
	{
	    row_sum += std::abs(H_values[j]);
	}

	for ( int j = 0; j < H_size; ++j )
	{
	    if ( row_sum > 0.0 )
	    {
		local_Q = std::abs(H_values[j]) / row_sum;
	    }
	    else
	    {
		local_Q = 0.0;
	    }
	    Q.InsertGlobalValues( H_indices[j], 1, &local_Q, &i );
	}
    }

    Q.FillComplete();
    Q.OptimizeStorage();
    return Q;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build the cumulative distribution function.
 */
Epetra_CrsMatrix AdjointMC::buildC()
{
    int N = d_Q.NumGlobalRows();
    int n_Q = d_Q.GlobalMaxNumEntries();
    Epetra_CrsMatrix C( Copy, d_Q.RowMap(), d_Q.GlobalMaxNumEntries() );
    double local_C = 0.0;
    std::vector<double> Q_values( n_Q );
    std::vector<int> Q_indices( n_Q );
    int size_Q = 0;
    for ( int i = 0; i < N; ++i )
    {
	d_Q.ExtractGlobalRowCopy( i, 
				  n_Q, 
				  size_Q, 
				  &Q_values[0], 
				  &Q_indices[0] );

	local_C = 0.0;
	for ( int j = 0; j < size_Q; ++j )
	{
	    local_C += Q_values[j];
	    C.InsertGlobalValues( i, 1, &local_C, &Q_indices[j] );
	}
    }

    C.FillComplete();
    C.OptimizeStorage();
    return C;
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Chimera_AdjointMC.cpp
//---------------------------------------------------------------------------//

