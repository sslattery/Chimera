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
 * \file Chimera_JacobiSplit_def.hpp
 * \author Stuart R. Slattery
 * \brief Jacobi split definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_JACOBISPLIT_DEF_HPP
#define Chimera_JACOBISPLIT_DEF_HPP

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_Map.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO>
JacobiSplit<Scalar,LO,GO>::JacobiSplit(
    const RCP_TpetraCrsMatrix& linear_operator )
{
    this->b_linear_operator = linear_operator;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Scalar, class LO, class GO>
JacobiSplit<Scalar,LO,GO>::~JacobiSplit()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Split the operator with Jacobi splitting.
 */
template<class Scalar, class LO, class GO>
void JacobiSplit<Scalar,LO,GO>::split()
{
    // Grab the row and column maps from the operator.
    Teuchos::RCP<const Tpetra::Map<LO,GO> > row_map =
	this->b_linear_operator->getRowMap();

    Teuchos::RCP<const Tpetra::Map<LO,GO> > col_map =
	this->b_linear_operator->getColMap();

    // Build the iteration matrix by extracting the diagonal and scaling by
    // its inverse.
    this->b_iteration_matrix = Teuchos::rcp( 
	new Tpetra::CrsMatrix<Scalar,LO,GO>( 
	    row_map, col_map,
	    this->b_linear_operator->getGlobalMaxNumRowEntries() ) );

    Teuchos::ArrayView<const LO> local_col_indices;
    Teuchos::ArrayView<const Scalar> local_values;

    Teuchos::Array<LO> diag_col_index(1,0);
    Teuchos::Array<Scalar> diag_zero(1,0.0);

    for ( LO row_index = row_map->getMinLocalIndex();
	  row_index < row_map->getMaxLocalIndex() + 1;
	  ++row_index )
    {
	this->b_linear_operator->getLocalRowView( 
	    row_index, local_col_indices, local_values );

	this->b_iteration_matrix->replaceLocalValues(
	    row_index, local_col_indices, local_values );

	diag_col_index[0] = 
	    col_map->getLocalElement( row_map->getGlobalElement( row_index ) );

	this->b_iteration_matrix->replaceLocalValues(
	    row_index, diag_col_index(), diag_zero() );
    }

    RCP_TpetraVector diagonal_inv = 
	Tpetra::createVector<Scalar,LO,GO>( row_map );
    this->b_linear_operator->getLocalDiagCopy( *diagonal_inv );
    diagonal_inv->reciprocal( *diagonal_inv );

    this->b_iteration_matrix->leftScale( *diagonal_inv );

    this->b_iteration_matrix->fillComplete();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Apply M^-1 to a vector (M^-1 x = y).
 */
template<class Scalar, class LO, class GO>
void JacobiSplit<Scalar,LO,GO>::applyInvM( const RCP_TpetraVector& x,
					   RCP_TpetraVector& y )
{
    Teuchos::RCP<const Tpetra::Map<LO,GO> > row_map =
	this->b_linear_operator->getRowMap();

    RCP_TpetraVector diagonal_inv = 
	Tpetra::createVector<Scalar,LO,GO>( row_map );
    this->b_linear_operator->getLocalDiagCopy( *diagonal_inv );
    diagonal_inv->reciprocal( *diagonal_inv );

    y->elementWiseMultiply( 0.0, *diagonal_inv, *x, 1.0 );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_JACOBISPLIT_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_JacobiSplit_def.hpp
//---------------------------------------------------------------------------//

