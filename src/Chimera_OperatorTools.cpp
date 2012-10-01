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
// \file Chimera_OperatorTools.cpp
// \author Stuart R. Slattery
// \brief OperatorTools definition.
//---------------------------------------------------------------------------//

#include "Chimera_OperatorTools.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>

#include <AnasaziTypes.hpp>
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"

namespace Chimera
{
namespace Solvers
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
OperatorTools::OperatorTools()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the spectral radius of an operator.
 */
double 
OperatorTools::spectralRadius( const Teuchos::RCP<Epetra_Operator>& matrix )
{
    Epetra_Map map = matrix->OperatorDomainMap();

    const int nev = 1;
    const int block_size = 1;
    const int num_blocks = 10;
    const int max_restarts = 100;
    const double tol = 1.0e-4;

    Teuchos::ParameterList krylovschur_params;
    krylovschur_params.set( "Which", "LM" );
    krylovschur_params.set( "Block Size", block_size );
    krylovschur_params.set( "Num Blocks", num_blocks );
    krylovschur_params.set( "Maximum Restarts", max_restarts );
    krylovschur_params.set( "Convergence Tolerance", tol );

    Teuchos::RCP<Epetra_MultiVector> vec 
	= Teuchos::rcp( new Epetra_MultiVector( map, block_size ) );
    vec->Random();

    Teuchos::RCP<Anasazi::BasicEigenproblem< double, 
					     Epetra_MultiVector, 
					     Epetra_Operator> > eigen_problem =
	Teuchos::rcp( 
	    new Anasazi::BasicEigenproblem< double, 
					    Epetra_MultiVector, 
					    Epetra_Operator >(matrix, vec) );
    eigen_problem->setNEV( nev );
    eigen_problem->setProblem();

    Anasazi::BlockKrylovSchurSolMgr< double, 
				     Epetra_MultiVector, 
				     Epetra_Operator> 
	solver_manager(eigen_problem, krylovschur_params);
    solver_manager.solve();

    Anasazi::Eigensolution<double,Epetra_MultiVector> sol = 
	eigen_problem->getSolution();
    std::vector<Anasazi::Value<double> > evals = sol.Evals;
    Teuchos::RCP<Epetra_MultiVector> evecs = sol.Evecs;

    double spectral_radius = pow( evals[0].realpart*evals[0].realpart +
				  evals[0].imagpart*evals[0].imagpart,
				  0.5 );

    return spectral_radius;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the stiffness ratio of an operator.
 */
double 
OperatorTools::stiffnessRatio( const Teuchos::RCP<Epetra_Operator>& matrix )
{
    return 0.0;
}

//---------------------------------------------------------------------------//

} // end namespace Solvers
} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Chimera_OperatorTools.cpp
//---------------------------------------------------------------------------//

