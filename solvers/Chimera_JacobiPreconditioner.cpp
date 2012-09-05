//---------------------------------------------------------------------------//
// \file JacobiPreconditioner.cpp
// \author Stuart R. Slattery
// \brief Jacobi preconditioner definition.
//---------------------------------------------------------------------------//

#include <vector>

#include "JacobiPreconditioner.hpp"

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

namespace HMCSA
{

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

/*!
 * \brief Destructor.
 */
JacobiPreconditioner::~JacobiPreconditioner()
{ /* ... */ }

/*!
 * \brief Do preconditioning on the operator and the right hand side.
 */
void JacobiPreconditioner::precondition()
{
    preconditionOperator();
    preconditionRHS();
}

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

} // end namespace HMCSA

//---------------------------------------------------------------------------//
// end JacobiPreconditioner.cpp
//---------------------------------------------------------------------------//

