//---------------------------------------------------------------------------//
// \file Chimera_MCSA.cpp
// \author Stuart R. Slattery
// \brief Monte Carlo Synthetic Acceleration solver definition.
//---------------------------------------------------------------------------//

#include "Chimera_MCSA.hpp"
#include "Chimera_AdjointMC.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

namespace Chimera
{
namespace Solvers
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor.
 */
MCSA::MCSA( Teuchos::RCP<Epetra_LinearProblem> &linear_problem )
    : d_linear_problem( linear_problem )
    , d_num_iters( 0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
MCSA::~MCSA()
{ /* ... */ }
 
//---------------------------------------------------------------------------//
/*!
 * \brief Solve.
 */
void MCSA::iterate( const int max_iters, const double tolerance,
		    const int num_histories, const double weight_cutoff )
{
    // Extract the linear problem.
    Epetra_CrsMatrix *A = 
	dynamic_cast<Epetra_CrsMatrix*>( d_linear_problem->GetMatrix() );
    Epetra_Vector *x = 
	dynamic_cast<Epetra_Vector*>( d_linear_problem->GetLHS() );
    const Epetra_Vector *b = 
	dynamic_cast<Epetra_Vector*>( d_linear_problem->GetRHS() );

    // Setup the residual Adjoint MC solver.
    Epetra_Map row_map = A->RowMap();
    Epetra_Vector delta_x( row_map );
    Epetra_Vector residual( row_map );
    Teuchos::RCP<Epetra_LinearProblem> residual_problem = Teuchos::rcp(
	new Epetra_LinearProblem( A, &delta_x, &residual ) );
    AdjointMC mc_solver = AdjointMC( residual_problem );

    // Iterate.
    Epetra_CrsMatrix H = mc_solver.getH();
    Epetra_Vector temp_vec( row_map );
    d_num_iters = 0;
    double residual_norm = 1.0;
    double b_norm;
    b->NormInf( &b_norm );
    double conv_crit = b_norm*tolerance;
    while ( residual_norm > conv_crit && d_num_iters < max_iters )
    {
	// Richardson iteration.
	H.Apply( *x, temp_vec );
	x->Update( 1.0, temp_vec, 1.0, *b, 0.0 );

	// Compute the residual.
	A->Apply( *x, temp_vec );
	residual.Update( 1.0, *b, -1.0, temp_vec, 0.0 );

	// Solve for delta_x.
	mc_solver.walk( num_histories, weight_cutoff );

	// Apply delta_x.
	x->Update( 1.0, delta_x, 1.0 );

	// Update convergence check.
	residual.NormInf( &residual_norm );
	++d_num_iters;
    }

    std::cout << "MCSA Iterations: " << d_num_iters << std::endl;
    std::cout << "Residual: " << residual_norm << std::endl;
}

//---------------------------------------------------------------------------//

} // end namespace Solvers
} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Chimera_MCSA.cpp
//---------------------------------------------------------------------------//

