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
// \file Chimera_MCSA.cpp
// \author Stuart R. Slattery
// \brief Monte Carlo Synthetic Acceleration solver definition.
//---------------------------------------------------------------------------//

#include "Chimera_MCSA.hpp"
#include "Chimera_AdjointMCDev.hpp"
#include "Chimera_DirectMC.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*! 
 * \brief Constructor.
 */
MCSA::MCSA( Teuchos::RCP<Epetra_LinearProblem> &linear_problem,
	    Teuchos::RCP<Teuchos::ParameterList> &plist )
    : d_linear_problem( linear_problem )
    , d_plist( plist )
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
void MCSA::iterate()
{
    // Get the solver parameters.
    int max_iters = d_plist->get<int>("MAX ITERS");
    double tolerance = d_plist->get<double>("TOLERANCE");

    // Extract the linear problem.
    Epetra_CrsMatrix *A = 
	dynamic_cast<Epetra_CrsMatrix*>( d_linear_problem->GetMatrix() );
    Epetra_Vector *x = 
	dynamic_cast<Epetra_Vector*>( d_linear_problem->GetLHS() );
    const Epetra_Vector *b = 
	dynamic_cast<Epetra_Vector*>( d_linear_problem->GetRHS() );

    // Setup the residual Adjoint MC solver.
    Epetra_Map row_map = A->RowMap();
    Epetra_Vector* delta_x = new Epetra_Vector( row_map );
    Epetra_Vector* residual = new Epetra_Vector( row_map );
    Teuchos::RCP<Epetra_LinearProblem> residual_problem = Teuchos::rcp(
	new Epetra_LinearProblem( A, delta_x, residual ) );
    AdjointMC mc_solver( residual_problem, d_plist );

    // Iterate.
    Epetra_CrsMatrix H = *mc_solver.getH();
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
	residual->Update( 1.0, *b, -1.0, temp_vec, 0.0 );

	// Solve for delta_x.
	delta_x->PutScalar( 0.0 );
	mc_solver.walk();

	// Apply delta_x.
	x->Update( 1.0, *delta_x, 1.0 );

	// Update convergence check.
	residual->NormInf( &residual_norm );
	++d_num_iters;

	std::cout << d_num_iters << " " << residual_norm << std::endl;
    }

    std::cout << "MCSA Iterations: " << d_num_iters << std::endl;
    std::cout << "Residual: " << residual_norm << std::endl;
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Chimera_MCSA.cpp
//---------------------------------------------------------------------------//

