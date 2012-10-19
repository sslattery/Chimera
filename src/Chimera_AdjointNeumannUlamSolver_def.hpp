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

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

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
{
    this->b_linear_problem = linear_problem;
    this->b_operator_split = operator_split;
    this->b_rng = rng;
    this->b_weight_cutoff = plist->get<Scalar>("WEIGHT CUTOFF");
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

};

} // end namespace Chimera

#endif // end Chimera_ADJOINTNEUMANNULAMSOLVER_HPP

//---------------------------------------------------------------------------//
// end Chimera_AdjointNeumannUlamSolver.hpp
//---------------------------------------------------------------------------//
