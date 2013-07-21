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

#include <algorithm>

#include "Chimera_Assertion.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Map.hpp>

#include <AnasaziTypes.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziBasicOutputManager.hpp>
#include <AnasaziTpetraAdapter.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Get a local component of an operator given a local row and column
 * index. 
 */
template<class Scalar, class LO, class GO>
Scalar OperatorTools::getMatrixComponentFromLocal( 
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& matrix,
    const LO local_row, const LO local_col )
{
    testPrecondition( matrix->getRowMap()->isNodeLocalElement( local_row ) );
    testPrecondition( matrix->getColMap()->isNodeLocalElement( local_col ) );

    Teuchos::ArrayView<const LO> local_indices;
    Teuchos::ArrayView<const Scalar> local_values;
    matrix->getLocalRowView( local_row, local_indices, local_values );

    typename Teuchos::ArrayView<const LO>::const_iterator local_idx_it =
	std::find( local_indices.begin(), local_indices.end(), local_col );

    if ( local_idx_it != local_indices.end() )
    {
	return local_values[ std::distance( local_indices.begin(),
					    local_idx_it ) ];
    }
    else
    {
	return 0.0;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get a local component of an operator given a global row and column
 * index. 
 */
template<class Scalar, class LO, class GO>
Scalar OperatorTools::getMatrixComponentFromGlobal( 
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& matrix,
    const GO global_row, const GO global_col )
{
    Teuchos::RCP<const Tpetra::Map<LO,GO> > row_map = matrix->getRowMap();
    Teuchos::RCP<const Tpetra::Map<LO,GO> > col_map = matrix->getColMap();

    LO local_row = row_map->getLocalElement( global_row );
    LO local_col = col_map->getLocalElement( global_col );

    testPrecondition( local_row != Teuchos::OrdinalTraits<LO>::invalid() );
    testPrecondition( local_col != Teuchos::OrdinalTraits<LO>::invalid() );

    Teuchos::ArrayView<const LO> local_indices;
    Teuchos::ArrayView<const Scalar> local_values;
    matrix->getLocalRowView( local_row, local_indices, local_values );

    typename Teuchos::ArrayView<const LO>::const_iterator local_idx_it =
	std::find( local_indices.begin(), local_indices.end(), local_col );

    if ( local_idx_it != local_indices.end() )
    {
	return local_values[ std::distance( local_indices.begin(),
					    local_idx_it ) ];
    }
    else
    {
	return 0.0;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Get the non-zero global column indices of a matrix that correspond
 *  to global row indices that are off-process. 
 */
template<class Scalar, class LO, class GO>
Teuchos::Array<GO> OperatorTools::getOffProcColumns(
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& matrix )
{
    Teuchos::RCP<const Tpetra::Map<LO,GO> > row_map = matrix->getRowMap();
    Teuchos::RCP<const Tpetra::Map<LO,GO> > col_map = matrix->getColMap();

    Teuchos::ArrayView<const GO> global_cols = col_map->getNodeElementList();
    typename Teuchos::ArrayView<const GO>::const_iterator global_col_it;

    Teuchos::Array<GO> off_proc_cols(0);

    for ( global_col_it = global_cols.begin();
	  global_col_it != global_cols.end();
	  ++global_col_it )
    {
	if ( !row_map->isNodeGlobalElement( *global_col_it ) )
	{
	    off_proc_cols.push_back( *global_col_it );
	}
    }

    return off_proc_cols;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the spectral radius of an operator.
 */
template<class Scalar,class LO, class GO>
void OperatorTools::spectralRadius( 
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& matrix,
    double& spec_rad,
    Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> >& eigenvector )
{
    typedef Tpetra::MultiVector<Scalar,LO,GO> MV;
    typedef Tpetra::Operator<Scalar,LO,GO> OP;

    const int nev = 1;
    const int block_size = 1;
    const int num_blocks = 10;
    const int max_restarts = 100;
    const Scalar tol = 1.0e-8;

    int verbosity = Anasazi::Errors + Anasazi::Warnings + 
                    Anasazi::FinalSummary + Anasazi::TimingDetails;

    Teuchos::ParameterList krylovschur_params;
    krylovschur_params.set( "Verbosity", verbosity );
    krylovschur_params.set( "Which", "LM" );
    krylovschur_params.set( "Block Size", block_size );
    krylovschur_params.set( "Num Blocks", num_blocks );
    krylovschur_params.set( "Maximum Restarts", max_restarts );
    krylovschur_params.set( "Convergence Tolerance", tol );

    Teuchos::RCP<MV> vec = 
	Teuchos::rcp( new MV(matrix->getRowMap(), block_size) );
    vec->randomize();

    Teuchos::RCP<Anasazi::BasicEigenproblem<Scalar,MV,OP> > eigen_problem =
	Teuchos::rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(
			  Teuchos::rcp_implicit_cast<OP>(matrix), vec ) );
    eigen_problem->setNEV( nev );
    eigen_problem->setProblem();

    Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP> solver_manager(
	eigen_problem, krylovschur_params );
    solver_manager.solve();

    Anasazi::Eigensolution<Scalar,MV > sol = eigen_problem->getSolution();
    std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
    Teuchos::RCP<MV> evecs = sol.Evecs;

    testPostcondition( sol.numVecs > 0 );

    spec_rad = std::pow( evals[0].realpart*evals[0].realpart +
                         evals[0].imagpart*evals[0].imagpart, 0.5 );

    eigenvector = evecs->getVectorNonConst(0);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the condition number of an operator.
 */
template<class Scalar,class LO, class GO>
Scalar OperatorTools::conditionNumber( 
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO> >& matrix )
{
    typedef Tpetra::MultiVector<Scalar,LO,GO> MV;
    typedef Tpetra::Operator<Scalar,LO,GO> OP;

    const int nev = 1;
    const int block_size = 1;
    const int num_blocks = 10;
    const int max_restarts = 100;
    const Scalar tol = 1.0e-4;

    Teuchos::ParameterList krylovschur_params;
    krylovschur_params.set( "Which", "SM" );
    krylovschur_params.set( "Block Size", block_size );
    krylovschur_params.set( "Num Blocks", num_blocks );
    krylovschur_params.set( "Maximum Restarts", max_restarts );
    krylovschur_params.set( "Convergence Tolerance", tol );

    Teuchos::RCP<MV> vec = 
	Teuchos::rcp( new MV(matrix->getRowMap(), block_size) );
    vec->randomize();

    Teuchos::RCP<Anasazi::BasicEigenproblem<Scalar,MV,OP> > eigen_problem =
	Teuchos::rcp( new Anasazi::BasicEigenproblem<Scalar,MV,OP>(
			  Teuchos::rcp_implicit_cast<OP>(matrix), vec ) );
    eigen_problem->setNEV( nev );
    eigen_problem->setProblem();

    Anasazi::BlockKrylovSchurSolMgr<Scalar,MV,OP> solver_manager(
	eigen_problem, krylovschur_params );
    solver_manager.solve();

    Anasazi::Eigensolution<Scalar,MV > sol = eigen_problem->getSolution();
    std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
    Teuchos::RCP<MV> evecs = sol.Evecs;

    testPostcondition( sol.numVecs > 0 );
    
    return std::pow( evals[0].realpart*evals[0].realpart +
		     evals[0].imagpart*evals[0].imagpart, 0.5 );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_OPERATORTOOLS_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_OperatorTools_def.cpp
//---------------------------------------------------------------------------//

