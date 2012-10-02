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
// \file Chimera_JacobiPreconditioner.cpp
// \author Stuart R. Slattery
// \brief Jacobi preconditioner definition.
//---------------------------------------------------------------------------//

#include <vector>

#include "Chimera_JacobiPreconditioner.hpp"

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
JacobiPreconditioner::JacobiPreconditioner( 
    Teuchos::RCP<Epetra_LinearProblem> &linear_problem )
    : d_linear_problem( linear_problem )
{
    // Get unconditioned operator.
    const Epetra_CrsMatrix *A = 
	dynamic_cast<Epetra_CrsMatrix*>( d_linear_problem->GetMatrix() );
    int N = A->NumGlobalRows();

    // Extract the diagonal.
    Epetra_Map map = A->RowMap();
    d_diagonal = Teuchos::rcp( new Epetra_Vector( map, false ) );
    A->ExtractDiagonalCopy( *d_diagonal );

    // Setup preconditioned system.
    d_M_inv_b = Teuchos::rcp(
	new Epetra_Vector( map, false ) );

    std::vector<int> entries_per_row( N, 1 );
    d_M_inv_A = Teuchos::rcp( 
	new Epetra_CrsMatrix( Copy, map, &entries_per_row[0] ) );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
JacobiPreconditioner::~JacobiPreconditioner()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Do preconditioning on the operator and the right hand side.
 */
void JacobiPreconditioner::precondition()
{
    preconditionOperator();
    preconditionRHS();
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do preconditioning on the operator.
 */
void JacobiPreconditioner::preconditionOperator()
{
    // Get unconditioned system.
    const Epetra_CrsMatrix *A = 
	dynamic_cast<Epetra_CrsMatrix*>( d_linear_problem->GetMatrix() );
    int N = A->NumGlobalRows();
    int n_A = A->GlobalMaxNumEntries();

    // Compute (M^-1 A) and (M^-1 b).
    double ma_val;
    int A_size;
    std::vector<double> A_values(n_A);
    std::vector<int> A_indices(n_A);
    for ( int i = 0; i < N; ++i )
    {
	A->ExtractGlobalRowCopy( i, n_A, A_size, &A_values[0], &A_indices[0] );
	for ( int j = 0; j < A_size; ++j )
	{
	    ma_val = A_values[j] / (*d_diagonal)[i];
	    d_M_inv_A->InsertGlobalValues( i, 1, &ma_val, &A_indices[j] );
	}
    }
    d_M_inv_A->FillComplete();
    d_M_inv_A->OptimizeStorage();

    // Modify the linear problem with the preconditioned operator.
    d_linear_problem->SetOperator( d_M_inv_A.getRawPtr() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Do preconditioning on the right hand side.
 */
void JacobiPreconditioner::preconditionRHS()
{
    // Get unconditioned system.
    const Epetra_Vector *b= 
	dynamic_cast<Epetra_Vector*>( d_linear_problem->GetRHS() );
    int N = b->GlobalLength();

    // Setup preconditioned system.
    std::vector<int> entries_per_row( N, 1 );

    // Compute (M^-1 b).
    for ( int i = 0; i < N; ++i )
    {
	(*d_M_inv_b)[i] = (*b)[i] / (*d_diagonal)[i];
    }

    // Modify the linear problem with the preconditioned rhs.
    d_linear_problem->SetRHS( d_M_inv_b.getRawPtr() );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Chimera_JacobiPreconditioner.cpp
//---------------------------------------------------------------------------//

