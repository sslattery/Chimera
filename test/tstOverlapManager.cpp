//---------------------------------------------------------------------------//
/*!
 * \file tstOverlapManager.cpp
 * \author Stuart R. Slattery
 * \brief OverlapManager tests.
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
#include <set>

#include <Chimera_LinearProblem.hpp>
#include <Chimera_OverlapManager.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_as.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

//---------------------------------------------------------------------------//
// Helper Functions.
//---------------------------------------------------------------------------//
// Build the linear Problem.
Teuchos::RCP<Chimera::LinearProblem<double,int,int> > buildLinearProblem()
{
    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    // Build the linear operator - this is a 2D Transient Diffusion operator.
    int N = 100;
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
    return Teuchos::rcp( new Chimera::LinearProblem<double,int,int>( A, X, B ) );
}

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
// Fully domain decomposed case.
TEUCHOS_UNIT_TEST( OverlapManager, no_overlap_test )
{
    // Set the overlap parameters.
    int num_overlap = 0;
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<int>("NUM OVERLAP", num_overlap);

    // Build the linear problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > linear_problem =
	buildLinearProblem();

    // Build the decomposition manager.
    Chimera::OverlapManager<double,int,int> overlap_manager( 
	linear_problem->getOperator(), linear_problem->getOperator(),
	linear_problem->getLHS(), plist );

    // Check the state parameters.
    TEST_ASSERT( overlap_manager.getNumOverlap() == num_overlap );

    // Get the base decomposition.
    Teuchos::RCP<const Tpetra::Map<int,int> > base_map = 
	linear_problem->getOperator()->getRowMap();
    Teuchos::ArrayView<const int> base_rows = 
	base_map->getNodeElementList();

    // Check that none of the base decomposition rows are in the overlap.
    for ( int i = 0; i < Teuchos::as<int>(base_rows.size()); ++i )
    {
	TEST_ASSERT( !overlap_manager.isOverlapGlobalElement( base_rows[i] ) );
    }

    // There is no overlap here so check that we have the same decomposition
    // as the input problem for the probability matrix.
    Teuchos::RCP<const Tpetra::Map<int,int> > probability_map = 
	overlap_manager.getOverlapProbabilityMatrix()->getRowMap();
    Teuchos::ArrayView<const int> probability_rows = 
	probability_map->getNodeElementList();

    TEST_ASSERT( base_rows.size() == probability_rows.size() );
    for ( int i = 0; i < Teuchos::as<int>(base_rows.size()); ++i )
    {
	TEST_ASSERT( base_rows[i] == probability_rows[i] );
    }

    // There is no overlap here so check that we have the same decomposition
    // as the input problem for the LHS.
    Teuchos::RCP<const Tpetra::Map<int,int> > lhs_map = 
	overlap_manager.getOverlapLHS()->getMap();
    Teuchos::ArrayView<const int> lhs_rows = 
	lhs_map->getNodeElementList();

    TEST_ASSERT( base_rows.size() == lhs_rows.size() );
    for ( int i = 0; i < Teuchos::as<int>(base_rows.size()); ++i )
    {
	TEST_ASSERT( base_rows[i] == lhs_rows[i] );
    }

    // Check that we have an additional state overlap for the iteration
    // matrix.
    Teuchos::RCP<const Tpetra::Map<int,int> > iteration_map = 
	overlap_manager.getOverlapIterationMatrix()->getRowMap();
    Teuchos::ArrayView<const int> iteration_rows = 
	iteration_map->getNodeElementList();

    std::set<int> ghosted_rows;
    for ( int i = 0; i < Teuchos::as<int>(iteration_rows.size()); ++i )
    {
	ghosted_rows.insert( iteration_rows[i] );
    }

    Teuchos::RCP<const Tpetra::Map<int,int> > iteration_col_map =
	overlap_manager.getOverlapIterationMatrix()->getColMap();
    Teuchos::ArrayView<const int> iteration_cols =
	iteration_col_map->getNodeElementList();

    for ( int i = 0; i < Teuchos::as<int>(iteration_cols.size()); ++i )
    {
	ghosted_rows.insert( iteration_cols[i] );
    }

    TEST_ASSERT( ghosted_rows.size() == 
		 Teuchos::as<std::size_t>(iteration_rows.size()) );

    Teuchos::ArrayView<const int>::const_iterator iteration_rows_it;
    std::set<int>::const_iterator ghosted_rows_it;
    for ( ghosted_rows_it = ghosted_rows.begin(), 
	iteration_rows_it = iteration_rows.begin();
	  ghosted_rows_it != ghosted_rows.end();
	  ++ghosted_rows_it, ++iteration_rows_it );
    {
	TEST_ASSERT( *ghosted_rows_it == *iteration_rows_it );
    }
}

//---------------------------------------------------------------------------//
// 1 state overlap.
TEUCHOS_UNIT_TEST( OverlapManager, overlap_1_test )
{
}

//---------------------------------------------------------------------------//
// 2 state overlap case.
TEUCHOS_UNIT_TEST( OverlapManager, overlap_2_test )
{
}

//---------------------------------------------------------------------------//
// end tstOverlapManager.cpp
//---------------------------------------------------------------------------//

