//---------------------------------------------------------------------------//
/*!
 * \file tstJacobiSplit.cpp
 * \author Stuart R. Slattery
 * \brief Boost random number generator unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <Chimera_LinearOperatorSplit.hpp>
#include <Chimera_JacobiSplit.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( JacobiSplit, jacobi_split_test )
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
    std::cout << "HERE" << std::endl;
    for ( int i = 0; i < local_num_rows; ++i )
    {
	gid = i + local_num_rows*comm_rank;
	global_columns[0] = gid-1;
	global_columns[1] = gid;
	global_columns[2] = gid+1;
	A->insertGlobalValues( gid, global_columns(), matrix_values() );
    }
    std::cout << "HERE" << std::endl;
    A->fillComplete();
    std::cout << "HERE" << std::endl;

    // Build the Jacobi split.
    Teuchos::RCP<Chimera::LinearOperatorSplit<double,int,int> > a_split =
	Teuchos::rcp( new Chimera::JacobiSplit<double,int,int>( A ) );
    a_split->split();

    // Check the operator.
    TEST_ASSERT( A == a_split->linearOperator() );

    // Check the iteration matrix.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int,int> > iteration_matrix =
	a_split->iterationMatrix();
    Teuchos::ArrayView<const double> row_values;
    Teuchos::ArrayView<const int> col_indices;
    for ( int i = 0; i < local_num_rows; ++i )
    {
	iteration_matrix->getLocalRowView( i, col_indices, row_values );

	TEST_ASSERT( row_values[0] == matrix_values[0] / matrix_values[1] );
	TEST_ASSERT( row_values[1] == 0.0 );
	TEST_ASSERT( row_values[2] == matrix_values[2] / matrix_values[1] );
    }

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
// end tstJacobiSplit.cpp
//---------------------------------------------------------------------------//

