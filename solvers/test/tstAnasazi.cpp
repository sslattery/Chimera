//---------------------------------------------------------------------------//
/*!
 * \file   mesh/test/tstAnasazi.cpp
 * \author Stuart Slattery
 * \brief  Anasazi class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <ostream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziEpetraAdapter.hpp"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( Anasazi, block_krylov_schur_test)
{
    int problem_size = 10;

    Epetra_SerialComm comm;
    Epetra_Map map( problem_size, 0, comm );

    // Build A.
    Teuchos::RCP<Epetra_CrsMatrix> A = 
	Teuchos::rcp( new Epetra_CrsMatrix( Copy, map, problem_size ) );
    double lower_diag = -1.0;
    double diag = 2.0;
    double upper_diag = -1.0;
    int global_row = 0;
    int lower_row = 0;
    int upper_row = 0;
    for ( int i = 0; i < problem_size; ++i )
    {
	global_row = A->GRID(i);
	lower_row = i-1;
	upper_row = i+1;
	if ( lower_row > -1 )
	{
	    A->InsertGlobalValues( global_row, 1, &lower_diag, &lower_row );
	}
	A->InsertGlobalValues( global_row, 1, &diag, &global_row );
	if ( upper_row < problem_size )
	{
	    A->InsertGlobalValues( global_row, 1, &upper_diag, &upper_row );
	}
    }
    A->FillComplete();

    // Block KrylovSchur setup.
    typedef Epetra_MultiVector MV;
    typedef Epetra_Operator OP;

    const int nev = 1;
    const int block_size = 2;
    const int num_blocks = 3;
    const int max_restarts = 100;
    const double tol = 1.0e-8;

    Teuchos::ParameterList solver_params;
    solver_params.set( "Which", "LM" );
    solver_params.set( "Block Size", block_size );
    solver_params.set( "Num Blocks", num_blocks );
    solver_params.set( "Maximum Restarts", max_restarts );
    solver_params.set( "Convergence Tolerance", tol );

    Teuchos::RCP<Epetra_MultiVector> ivec 
	= Teuchos::rcp( new Epetra_MultiVector(map, block_size) );
    ivec->Random();

    // Create the eigenproblem.
    Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
	Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );

    // Set the number of eigenvalues requested.
    MyProblem->setNEV( nev );

    // Finalize eigenproblem.
    bool boolret = MyProblem->setProblem();
    TEST_ASSERT( boolret );

    // Create the solver manager
    Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMan(MyProblem, 
							     solver_params);

    // Solve the problem
    MySolverMan.solve();

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
    std::vector<Anasazi::Value<double> > evals = sol.Evals;
    Teuchos::RCP<MV> evecs = sol.Evecs;

    std::cout << "NUM EIGV " << evals.size() << std::endl;
    for( int i = 0; i < (int) evals.size(); ++i )
    {
	std::cout << "EIGVAL " << i << " " 
		  << evals[i].realpart << " "
		  << evals[i].imagpart << std::endl;
    }
}

//---------------------------------------------------------------------------//
//                        end of tstAnasazi.cpp
//---------------------------------------------------------------------------//
