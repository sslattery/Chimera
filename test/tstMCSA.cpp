//---------------------------------------------------------------------------//
/*!
 * \file tstMCSA.cpp
 * \author Stuart R. Slattery
 * \brief MCSA tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <string>
#include <cassert>

#include <Chimera_LinearSolver.hpp>
#include <Chimera_LinearProblem.hpp>
#include <Chimera_LinearSolverFactory.hpp>
#include <Chimera_MCSA.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( MCSA, mcsa_test )
{
    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Build the linear operator - this is a 2D Transient Diffusion operator.
    int N = 10;
    int problem_size = N*N;
    double dx = 0.01;
    double dy = 0.01;
    double dt = 0.05;
    double alpha = 0.001;
    Teuchos::RCP<const Tpetra::Map<int> > row_map = 
	Tpetra::createUniformContigMap<int,int>( problem_size, comm );
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int> > A = 
	Tpetra::createCrsMatrix<double,int,int>( row_map, 3 );

    Teuchos::Array<double> diag( 1, 1.0 + 2*dt*alpha*( 1/(dx*dx) + 1/(dy*dy) ) );
    Teuchos::Array<double> i_minus( 1, -dt*alpha/(dx*dx) );
    Teuchos::Array<double> i_plus( 1, -dt*alpha/(dx*dx) );
    Teuchos::Array<double> j_minus( 1, -dt*alpha/(dy*dy) );
    Teuchos::Array<double> j_plus( 1, -dt*alpha/(dy*dy) );

    Teuchos::Array<int> idx(1);
    Teuchos::Array<int> idx_iminus(1);
    Teuchos::Array<int> idx_iplus(1);
    Teuchos::Array<int> idx_jminus(1);
    Teuchos::Array<int> idx_jplus(1);
    Teuchos::Array<double> one(1,1.0);

    if ( comm->getRank() == 0 )
    {
	// Min X boundary Dirichlet.
	for ( int j = 1; j < N-1; ++j )
	{
	    int i = 0;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Max X boundary Dirichlet.
	for ( int j = 1; j < N-1; ++j )
	{
	    int i = N-1;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Min Y boundary Dirichlet.
	for ( int i = 0; i < N; ++i )
	{
	    int j = 0;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Max Y boundary Dirichlet.
	for ( int i = 0; i < N; ++i )
	{
	    int j = N-1;
	    idx[0] = i + j*N;
	    A->insertGlobalValues( idx[0], idx(), one() );
	}

	// Central grid points.
	for ( int i = 1; i < N-1; ++i )
	{
	    for ( int j = 1; j < N-1; ++j )
	    {
		idx[0] = i + j*N;
		idx_iminus[0] = (i-1) + j*N;
		idx_iplus[0] = (i+1) + j*N;
		idx_jminus[0] = i + (j-1)*N;
		idx_jplus[0] = i + (j+1)*N;
	    
		A->insertGlobalValues( idx[0], idx_jminus(), j_minus() );
		A->insertGlobalValues( idx[0], idx_iminus(), i_minus() );
		A->insertGlobalValues( idx[0], idx(),        diag()    );
		A->insertGlobalValues( idx[0], idx_iplus(),  i_plus()  );
		A->insertGlobalValues( idx[0], idx_jplus(),  j_plus()  );
	    }
	}
    }
    comm->barrier();

    A->fillComplete();

    // Build the solution vector.
    double X_val = 0.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > X = 
	Tpetra::createVector<double,int>( row_map );
    X->putScalar( X_val );

    // Build the right hand side.
    double B_val = 0.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > B = 
	Tpetra::createVector<double,int>( row_map );
    B->putScalar( B_val );

    // Set dirichlet boundary conditions.
    int row;
    // left
    for ( int j = 1; j < N-1; ++j )
    {
	int i = 0;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 5.0 );
	}
    }
    // right
    for ( int j = 1; j < N-1; ++j )
    {
	int i = N-1;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 5.0 );
	}
    }
    // bottom
    for ( int i = 0; i < N; ++i )
    {
	int j = 0;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 0.0 );
	}
    }
    // top
    for ( int i = 0; i < N; ++i )
    {
	int j = N-1;
	row = i + j*N;
	if ( row_map->isNodeGlobalElement( row ) )
	{
	    B->replaceGlobalValue( row, 0.0 );
	}
    }

    // Build the linear problem.
    comm->barrier();
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > linear_problem = 
	Teuchos::rcp( new Chimera::LinearProblem<double,int,int>( A, X, B ) );

    // Set the solver parameters.
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<std::string>( "SOLVER TYPE",         "MCSA"    );
    plist->set<std::string>( "RNG TYPE",            "MT19937" );
    plist->set<double>(      "TOLERANCE",           1.0e-8    );
    plist->set<int>(         "MAX ITERS",           100       );
    plist->set<std::string>( "SPLIT TYPE",          "JACOBI"  );
    plist->set<std::string>( "MC TYPE",             "ADJOINT" );
    plist->set<double>(      "WEIGHT CUTOFF",       1.0e-4    );
    plist->set<int>(         "HISTORIES PER STAGE", 1000      );

    Teuchos::RCP<Chimera::LinearSolver<double,int,int> > solver =
	Chimera::LinearSolverFactory::create( plist, linear_problem );

    // Solve.
    solver->iterate();

    // Check the solution vector.
    Teuchos::ArrayRCP<const double> X_view = 
	linear_problem->getLHS()->get1dView();
    std::cout << X_view() << std::endl;
}

//---------------------------------------------------------------------------//
// end tstMCSA.cpp
//---------------------------------------------------------------------------//

