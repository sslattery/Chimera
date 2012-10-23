//---------------------------------------------------------------------------//
/*!
 * \file tstAdjointNeumannUlamSolver.cpp
 * \author Stuart R. Slattery
 * \brief Stationary solver tests.
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
#include <Chimera_LinearOperatorSplit.hpp>
#include <Chimera_LinearOperatorSplitFactory.hpp>
#include <Chimera_NeumannUlamSolver.hpp>
#include <Chimera_AdjointNeumannUlamSolver.hpp>
#include <Chimera_BoostRNG.hpp>
#include <Chimera_OperatorTools.hpp>
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
TEUCHOS_UNIT_TEST( AdjointNeumannUlamSolver, adjoint_neumannulam_test )
{
    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Build the linear operator - this is a 2D Transient Diffusion operator.
    int N = 10;
    int problem_size = N*N;
    double dx = 0.01;
    double dy = 0.01;
    double dt = 0.005;
    double alpha = 0.01;
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
	B->replaceGlobalValue( row, 5.0 );
    }
    // right
    for ( int j = 1; j < N-1; ++j )
    {
	int i = N-1;
	row = i + j*N;
	B->replaceGlobalValue( row, 5.0 );
    }
    // bottom
    for ( int i = 0; i < N; ++i )
    {
	int j = 0;
	row = i + j*N;
	B->replaceGlobalValue( row, 0.0 );
    }
    // top
    for ( int i = 0; i < N; ++i )
    {
	int j = N-1;
	row = i + j*N;
	B->replaceGlobalValue( row, 0.0 );
    }

    // Build the linear problem.
    comm->barrier();
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > linear_problem = 
	Teuchos::rcp( new Chimera::LinearProblem<double,int,int>( A, X, B ) );

    // Build the random number generator.
    Teuchos::RCP<boost::mt19937> rng = 
	Chimera::RNGTraits<boost::mt19937>::create();

    // Build the Adjoint solver.
    std::string split_type = "JACOBI";
    double weight_cutoff = 1.0e-4;
    int histories_per_stage = 1000;
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<std::string>("SPLIT TYPE", split_type);
    plist->set<double>("WEIGHT CUTOFF", weight_cutoff);
    plist->set<int>("HISTORIES PER STAGE", histories_per_stage);

    Teuchos::RCP<Chimera::LinearOperatorSplit<double,int,int> > lin_op_split =
	Chimera::LinearOperatorSplitFactory::create( plist, A );
    lin_op_split->split();

    Teuchos::RCP<Chimera::NeumannUlamSolver<double,int,int,boost::mt19937> > solver =
	Teuchos::rcp( 
	    new Chimera::AdjointNeumannUlamSolver<double,int,int,boost::mt19937>(
		linear_problem, lin_op_split, rng, plist ) );

    // Check the solver parameters.
    TEST_ASSERT( solver->linearProblem() == linear_problem );
    TEST_ASSERT( solver->linearOperatorSplit() == lin_op_split );
    TEST_ASSERT( solver->weightCutoff() == weight_cutoff );
    TEST_ASSERT( solver->historiesPerStage() == histories_per_stage );
    std::cout << "SPEC RAD " << Chimera::OperatorTools::spectralRadius(
	lin_op_split->iterationMatrix() ) << std::endl;
    // Solve.
    solver->walk();

    // Check the solution vector.
    Teuchos::ArrayRCP<const double> X_view = 
	linear_problem->getLHS()->get1dView();
    std::cout << X_view() << std::endl;
}

//---------------------------------------------------------------------------//
// end tstAdjointNeumannUlamSolver.cpp
//---------------------------------------------------------------------------//

