//---------------------------------------------------------------------------//
/*!
 * \file tstLinearOperatorSplitFactory.cpp
 * \author Stuart R. Slattery
 * \brief Linear operator split factory tests.
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

#include <Chimera_LinearOperatorSplit.hpp>
#include <Chimera_LinearOperatorSplitFactory.hpp>

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
TEUCHOS_UNIT_TEST( LinearOperatorSplitFactory, jacobi_split_test )
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
    matrix_values[2] = 3.0;
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
    edge_matrix_values[0] = 2.0;
    edge_matrix_values[1] = 3.0;
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

    edge_matrix_values[0] = 1.0;
    edge_matrix_values[1] = 2.0;
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

    // Build the operator splitting.
    Teuchos::RCP<Teuchos::ParameterList> plist = 
	Teuchos::rcp( new Teuchos::ParameterList() );
    plist->set<std::string>("SPLIT TYPE","JACOBI");

    Teuchos::RCP<Chimera::LinearOperatorSplit<double,int,int> > a_split =
	Chimera::LinearOperatorSplitFactory::create( plist, A );

    // Check the operator.
    TEST_ASSERT( A == a_split->linearOperator() );

    // Check the iteration matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int> > iteration_matrix =
	a_split->iterationMatrix();
    Teuchos::ArrayView<const double> row_values;
    Teuchos::ArrayView<const int> col_indices;
    for ( int i = 1; i < local_num_rows-1; ++i )
    {
	iteration_matrix->getLocalRowView( i, col_indices, row_values );
	TEST_ASSERT( row_values[0] == -matrix_values[0] / matrix_values[1] );
	TEST_ASSERT( row_values[1] == 0.0 );
	TEST_ASSERT( row_values[2] == -matrix_values[2] / matrix_values[1] );
    }

    iteration_matrix->getLocalRowView( 0, col_indices, row_values );
    if ( comm_rank == 0 )
    {
	TEST_ASSERT( row_values[0] == 0.0 );
	TEST_ASSERT( row_values[1] == -matrix_values[2] / matrix_values[1] );
    }
    else
    {
	TEST_ASSERT( row_values[0] == -matrix_values[0] / matrix_values[1] );
	TEST_ASSERT( row_values[1] == 0.0 );
	TEST_ASSERT( row_values[2] == -matrix_values[2] / matrix_values[1] );
    }
    comm->barrier();

    iteration_matrix->getLocalRowView( local_num_rows-1, col_indices, row_values );
    if ( comm_rank == comm_size-1 )
    {
	TEST_ASSERT( row_values[0] == -matrix_values[0] / matrix_values[1] );
	TEST_ASSERT( row_values[1] == 0.0 );
    }
    else
    {
	TEST_ASSERT( row_values[0] == -matrix_values[0] / matrix_values[1] );
	TEST_ASSERT( row_values[1] == 0.0 );
	TEST_ASSERT( row_values[2] == -matrix_values[2] / matrix_values[1] );
    }
    comm->barrier();

    // Check application to a vector where the product is another vector.
    double X_val = 5.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > X = 
	Tpetra::createVector<double,int>( row_map );
    X->putScalar( X_val );

    Teuchos::RCP<Tpetra::Vector<double,int> > Y = 
	Tpetra::createVector<double,int>( row_map );

    a_split->applyInvM( X, Y );
    Teuchos::ArrayRCP<const double> Y_view = Y->get1dView();
    for ( int i = 0; i < local_num_rows; ++i )
    {
	TEST_ASSERT( Y_view[i] == X_val / matrix_values[1] );
    }

    // Check the application to a vector where the product is itself.
    a_split->applyInvM( X, X );
    Teuchos::ArrayRCP<const double> X_view = X->get1dView();
    for ( int i = 0; i < local_num_rows; ++i )
    {
	TEST_ASSERT( X_view[i] == X_val / matrix_values[1] );
    }
}

//---------------------------------------------------------------------------//
// end tstLinearOperatorSplitFactory.cpp
//---------------------------------------------------------------------------//

