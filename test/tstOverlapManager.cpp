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
    // Set the decomposition parameters.
    double overlap_frac = 0.0;
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<double>("OVERLAP FRACTION", overlap_frac);

    // Build the linear problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > linear_problem =
	buildLinearProblem();

    // Build the decomposition manager.
    Chimera::OverlapManager<double,int,int> decomp_manager( 
	linear_problem, plist );

    // Check the state parameters.
    TEST_ASSERT( decomp_manager.getOverlapFraction() == overlap_frac );

    // Get the decomposed problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > decomp_problem =
	decomp_manager.getDecomposedProblem();

    // There is no overlap here so check that we have the same decomposition
    // as the input problem for the operator.
    Teuchos::RCP<const Tpetra::Map<int,int> > base_op_map = 
	linear_problem->getOperator()->getRowMap();
    Teuchos::ArrayView<const int> base_op_rows = 
	base_op_map->getNodeElementList();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_op_map = 
	decomp_problem->getOperator()->getRowMap();
    Teuchos::ArrayView<const int> decomp_op_rows = 
	decomp_op_map->getNodeElementList();

    TEST_ASSERT( base_op_rows.size() == decomp_op_rows.size() );
    for ( int i = 0; i < base_op_rows.size(); ++i )
    {
	TEST_ASSERT( base_op_rows[i] == decomp_op_rows[i] );
    }

    // There is no overlap here so check that we have the same decomposition
    // as the input problem for the LHS.
    Teuchos::RCP<const Tpetra::Map<int,int> > base_lhs_map = 
	linear_problem->getLHS()->getMap();
    Teuchos::ArrayView<const int> base_lhs_rows = 
	base_lhs_map->getNodeElementList();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_lhs_map = 
	decomp_problem->getLHS()->getMap();
    Teuchos::ArrayView<const int> decomp_lhs_rows = 
	decomp_lhs_map->getNodeElementList();

    TEST_ASSERT( base_lhs_rows.size() == decomp_lhs_rows.size() );
    for ( int i = 0; i < base_lhs_rows.size(); ++i )
    {
	TEST_ASSERT( base_lhs_rows[i] == decomp_lhs_rows[i] );
    }

    // The RHS should be the same as the input.
    Teuchos::RCP<const Tpetra::Map<int,int> > base_lhs_map = 
	linear_problem->getRHS()->getMap();
    Teuchos::ArrayView<const int> base_lhs_rows = 
	base_lhs_map->getNodeElementList();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_lhs_map = 
	decomp_problem->getRHS()->getMap();
    Teuchos::ArrayView<const int> decomp_lhs_rows = 
	decomp_lhs_map->getNodeElementList();

    TEST_ASSERT( base_lhs_rows.size() == decomp_lhs_rows.size() );
    for ( int i = 0; i < base_lhs_rows.size(); ++i )
    {
	TEST_ASSERT( base_lhs_rows[i] == decomp_lhs_rows[i] );
    }
}

//---------------------------------------------------------------------------//
// Fully replicated case.
TEUCHOS_UNIT_TEST( OverlapManager, full_overlap_test )
{
    // Set the decomposition parameters.
    double overlap_frac = 1.0;
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<double>("OVERLAP FRACTION", overlap_frac);

    // Build the linear problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > linear_problem =
	buildLinearProblem();

    // Build the decomposition manager.
    Chimera::OverlapManager<double,int,int> decomp_manager( 
	linear_problem, plist );

    // Check the state parameters.
    TEST_ASSERT( decomp_manager.getOverlapFraction() == overlap_frac );

    // Get the decomposed problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > decomp_problem =
	decomp_manager.getDecomposedProblem();

    // There is no overlap here so check that we have the full operator on all
    // procs.
    int base_op_num_rows = linear_problem->getOperator()->getGlobalNumRows();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_op_map = 
	decomp_problem->getOperator()->getRowMap();
    Teuchos::ArrayView<const int> decomp_op_rows = 
	decomp_op_map->getNodeElementList();

    TEST_ASSERT( base_op_num_rows == decomp_op_rows.size() );
    for ( int i = 0; i < base_op_num_rows; ++i )
    {
	TEST_ASSERT( i == decomp_op_rows[i] );
    }

    // There is full overlap here so check that we have the full LHS on all
    // procs. 
    int base_lhs_num_rows = linear_problem->getLHS()->getGlobalLength();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_lhs_map = 
	decomp_problem->getLHS()->getMap();
    Teuchos::ArrayView<const int> decomp_lhs_rows = 
	decomp_lhs_map->getNodeElementList();

    TEST_ASSERT( base_lhs_num_rows == decomp_lhs_rows.size() );
    for ( int i = 0; i < base_lhs_num_rows(); ++i )
    {
	TEST_ASSERT( i == decomp_lhs_rows[i] );
    }

    // The RHS should be the same as the input.
    Teuchos::RCP<const Tpetra::Map<int,int> > base_lhs_map = 
	linear_problem->getRHS()->getMap();
    Teuchos::ArrayView<const int> base_lhs_rows = 
	base_lhs_map->getNodeElementList();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_lhs_map = 
	decomp_problem->getRHS()->getMap();
    Teuchos::ArrayView<const int> decomp_lhs_rows = 
	decomp_lhs_map->getNodeElementList();

    TEST_ASSERT( base_lhs_rows.size() == decomp_lhs_rows.size() );
    for ( int i = 0; i < base_lhs_rows.size(); ++i )
    {
	TEST_ASSERT( base_lhs_rows[i] == decomp_lhs_rows[i] );
    }
}

//---------------------------------------------------------------------------//
// Half overlap decomposed case.
TEUCHOS_UNIT_TEST( OverlapManager, half_overlap_test )
{
    // Set the decomposition parameters.
    double overlap_frac = 0.5;
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<double>("OVERLAP FRACTION", overlap_frac);

    // Build the linear problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > linear_problem =
	buildLinearProblem();

    // Build the decomposition manager.
    Chimera::OverlapManager<double,int,int> decomp_manager( 
	linear_problem, plist );

    // Check the state parameters.
    TEST_ASSERT( decomp_manager.getOverlapFraction() == overlap_frac );

    // Get the decomposed problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > decomp_problem =
	decomp_manager.getDecomposedProblem();

    // There is half overlap here so check that we have the right operator
    // rows on proc.
    int base_op_num_rows = linear_problem->getOperator()->getGlobalNumRows();
    Teuchos::RCP<const Tpetra::Map<int,int> > base_op_map = 
	linear_problem->getOperator()->getRowMap();
    Teuchos::ArrayView<const int> base_op_rows = 
	base_op_map->getNodeElementList();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_op_map = 
	decomp_problem->getOperator()->getRowMap();
    Teuchos::ArrayView<const int> decomp_op_rows = 
	decomp_op_map->getNodeElementList();

    int num_overlap_rows = 
	std::floor(overlap_frac*base_op_rows.size()/base_op_num_rows);

    TEST_ASSERT( num_overlap_rows == decomp_op_rows.size() );
    for ( int i = 0; i < base_op_rows.size(); ++i )
    {
	TEST_ASSERT( base_op_rows[i] == decomp_op_rows[i] );
    }

    // There is half overlap here so check that we have the right LHS rows on
    // proc.
    int base_lhs_num_rows = linear_problem->getLHS()->getGlobalLength();
    Teuchos::RCP<const Tpetra::Map<int,int> > base_lhs_map = 
	linear_problem->getLHS()->getMap();
    Teuchos::ArrayView<const int> base_lhs_rows = 
	base_lhs_map->getNodeElementList();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_lhs_map = 
	decomp_problem->getLHS()->getMap();
    Teuchos::ArrayView<const int> decomp_lhs_rows = 
	decomp_lhs_map->getNodeElementList();

    TEST_ASSERT( num_overlap_rows == decomp_lhs_rows.size() );
    for ( int i = 0; i < base_lhs_rows.size(); ++i )
    {
	TEST_ASSERT( base_lhs_rows[i] == decomp_lhs_rows[i] );
    }

    // The RHS should be the same as the input.
    Teuchos::RCP<const Tpetra::Map<int,int> > base_lhs_map = 
	linear_problem->getRHS()->getMap();
    Teuchos::ArrayView<const int> base_lhs_rows = 
	base_lhs_map->getNodeElementList();

    Teuchos::RCP<const Tpetra::Map<int,int> > decomp_lhs_map = 
	decomp_problem->getRHS()->getMap();
    Teuchos::ArrayView<const int> decomp_lhs_rows = 
	decomp_lhs_map->getNodeElementList();

    TEST_ASSERT( base_lhs_rows.size() == decomp_lhs_rows.size() );
    for ( int i = 0; i < base_lhs_rows.size(); ++i )
    {
	TEST_ASSERT( base_lhs_rows[i] == decomp_lhs_rows[i] );
    }
}

//---------------------------------------------------------------------------//
// end tstOverlapManager.cpp
//---------------------------------------------------------------------------//

