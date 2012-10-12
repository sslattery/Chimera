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
// \file Chimera_OperatorTools_def.hpp
// \author Stuart R. Slattery
// \brief OperatorTools definition.
//---------------------------------------------------------------------------//

#ifndef Chimera_OPERATORTOOLS_DEF_HPP
#define Chimera_OPERATORTOOLS_DEF_HPP

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_MultiVector.hpp>

#include <AnasaziTypes.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziEpetraAdapter.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Compute the spectral radius of an operator.
 */
template<class Scalar,class LO, class GO>
Scalar OperatorTools::spectralRadius( 
    const Teuchos::RCP<Tpetra::Operator<Scalar,LO,GO> >& matrix )
{
    typedef Tpetra::MultiVector<Scalar,LO,GO> MV;
    typedef Tpetra::Operator<Scalar,LO,GO> OP;

    const int nev = 1;
    const int block_size = 1;
    const int num_blocks = 10;
    const int max_restarts = 100;
    const Scalar tol = 1.0e-4;

    Teuchos::ParameterList krylovschur_params;
    krylovschur_params.set( "Which", "LM" );
    krylovschur_params.set( "Block Size", block_size );
    krylovschur_params.set( "Num Blocks", num_blocks );
    krylovschur_params.set( "Maximum Restarts", max_restarts );
    krylovschur_params.set( "Convergence Tolerance", tol );

    Teuchos::RCP<MV > vec = Teuchos::rcp( 
	new MV(matrix->getRowMap, block_size) );
    vec->random();

    Teuchos::RCP<Anasazi::BasicEigenproblem<Scalar,MV,OP> > eigen_problem =
	Teuchos::rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(
			  matrix, vec ) );
    eigen_problem->setNEV( nev );
    eigen_problem->setProblem();

    Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP> solver_manager(
	eigen_problem, krylovschur_params );
    solver_manager.solve();

    Anasazi::Eigensolution<Scalar,MV > sol = eigen_problem->getSolution();
    std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;

    return pow( evals[0].realpart*evals[0].realpart +
		evals[0].imagpart*evals[0].imagpart, 0.5 );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_OPERATORTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_OperatorTools_def.cpp
//---------------------------------------------------------------------------//

