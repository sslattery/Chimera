//---------------------------------------------------------------------------//
/*!
 * \file tstStationarySolver.cpp
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
#include <Chimera_StationarySolver.hpp>
#include <Chimera_LinearProblem.hpp>

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
TEUCHOS_UNIT_TEST( StationarySolver, jacobi_solve_test )
{
    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Setup linear operator distribution.
    int local_num_rows = 10;
    int global_num_rows = local_num_rows*comm_size;
    Teuchos::RCP<const Tpetra::Map<int> > row_map = 
	Tpetra::createUniformContigMap<int,int>( global_num_rows, comm );

    // Build the linear operator.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int> > A = 
	Tpetra::createCrsMatrix<double,int,int>( row_map, 3 );

    Teuchos::Array<int> global_columns(3);
    Teuchos::Array<double> matrix_values(3);
    matrix_values[0] = 1.0;
    matrix_values[1] = 2.0;
    matrix_values[2] = 4.0;
    int gid = 0;
    for ( int i = 1; i < local_num_rows-1; ++i )
    {
	gid = i + local_num_rows*comm_rank;
	global_columns[0] = gid-1;
	global_columns[1] = gid;
	global_columns[2] = gid+1;
	A->insertGlobalValues( gid, global_columns(), matrix_values() );
    }

    Teuchos::Array<int> edge_global_columns(2);
    Teuchos::Array<double> edge_matrix_values(2);
    edge_matrix_values[0] = matrix_values[1];
    edge_matrix_values[1] = matrix_values[2];
    gid = local_num_rows*comm_rank;
    edge_global_columns[0] = gid;
    edge_global_columns[1] = gid+1;
    global_columns[0] = gid-1;
    global_columns[1] = gid;
    global_columns[2] = gid+1;
    if ( comm_rank == 0 )
    {
	A->insertGlobalValues( gid, edge_global_columns(), edge_matrix_values() );
    }
    else
    {
	A->insertGlobalValues( gid, global_columns(), matrix_values() );
    }
    comm->barrier();

    edge_matrix_values[0] = matrix_values[0];
    edge_matrix_values[1] = matrix_values[1];
    gid = local_num_rows - 1 + local_num_rows*comm_rank;
    edge_global_columns[0] = gid-1;
    edge_global_columns[1] = gid;
    global_columns[0] = gid-1;
    global_columns[1] = gid;
    global_columns[2] = gid+1;
    if ( comm_rank == comm_size-1 )
    {
	A->insertGlobalValues( gid, edge_global_columns(), edge_matrix_values() );
    }
    else
    {
	A->insertGlobalValues( gid, global_columns(), matrix_values() );
    }
    comm->barrier();

    A->fillComplete();

    // Build the solution vector.
    double X_val = 0.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > X = 
	Tpetra::createVector<double,int>( row_map );
    X->putScalar( X_val );

    // Build the right hand side.
    double B_val = 5.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > B = 
	Tpetra::createVector<double,int>( row_map );
    B->putScalar( B_val );

    // Build the linear problem.
    comm->barrier();
    Teuchos::RCP<Chimera::LinearProblem<double,int,int> > linear_problem = 
	Teuchos::rcp( new Chimera::LinearProblem<double,int,int>( A, X, B ) );

    // Build the stationary solver.
    std::string split_type = "JACOBI";
    double tolerance = 1.0e-8;
    int max_iters = 2;
    Teuchos::RCP<Teuchos::ParameterList> plist =
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<std::string>("SPLIT TYPE", split_type);
    plist->set<double>("TOLERANCE", tolerance);
    plist->set<int>("MAX ITERS", max_iters);
    Teuchos::RCP<Chimera::LinearSolver<double,int,int> > solver =
	Teuchos::rcp( new Chimera::StationarySolver<double,int,int>(
			  linear_problem, plist ) );

    // Check the solver parameters.
    TEST_ASSERT( solver->linearProblem() == linear_problem );
    TEST_ASSERT( solver->tolerance() == tolerance );
    TEST_ASSERT( solver->maxNumIters() == max_iters );

    // Solve (only 2 iterations).
    solver->iterate();

    // Check the solution vector.
    Teuchos::ArrayRCP<const double> X_view = linear_problem->getLHS()->get1dView();
    for ( int i = 1; i < local_num_rows-1; ++i )
    {
	TEST_ASSERT( X_view[i] == 
		     (1-(matrix_values[0]+matrix_values[2])/
		      matrix_values[1])*B_val/matrix_values[1] );
    }

    if ( comm_rank == 0 )
    {
	TEST_ASSERT( X_view[0] == 
		     (1-(matrix_values[2])/
		      matrix_values[1])*B_val/matrix_values[1] );   
    }
    else
    {
	TEST_ASSERT( X_view[0] == 
		     (1-(matrix_values[0]+matrix_values[2])/
		      matrix_values[1])*B_val/matrix_values[1] );   
    }
    comm->barrier();

    if ( comm_rank == comm_size-1 )
    {
	TEST_ASSERT( X_view[local_num_rows-1] == 
		     (1-(matrix_values[0])/
		      matrix_values[1])*B_val/matrix_values[1] );   
    }
    else
    {
	TEST_ASSERT( X_view[local_num_rows-1] == 
		     (1-(matrix_values[0]+matrix_values[2])/
		      matrix_values[1])*B_val/matrix_values[1] );   
    }
    comm->barrier();
}

//---------------------------------------------------------------------------//
// end tstStationarySolver.cpp
//---------------------------------------------------------------------------//

