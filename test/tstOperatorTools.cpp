//---------------------------------------------------------------------------//
/*!
 * \file tstOperatorTools.cpp
 * \author Stuart R. Slattery
 * \brief Operator tools tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <Chimera_OperatorTools.hpp>
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
#include <Tpetra_RowMatrix.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( OperatorTools, matrix_component_test )
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
    matrix_values[0] = -1.0;
    matrix_values[1] = 2.0;
    matrix_values[2] = -3.0;
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

    // Check that we can get the right components of the matrix.
    double matrix_component = 0.0;
    int global_i = 0;
    for ( int i = 0; i < local_num_rows; ++i )
    {
	for ( int j = 0; j < local_num_rows; ++j )
	{
	    global_i = comm_rank*local_num_rows + i;

	    matrix_component = 
		Chimera::OperatorTools::getMatrixComponentFromGlobal( 
		    A, global_i, j );

	    if ( global_i == j )
	    {
		TEST_ASSERT( matrix_values[1] == matrix_component );
	    }
	    else if ( global_i-1 == j )
	    {
		TEST_ASSERT( matrix_values[0] == matrix_component );
	    }
	    else if ( global_i+1 == j )
	    {
		TEST_ASSERT( matrix_values[2] == matrix_component );
	    }
	    else
	    {
		TEST_ASSERT( 0.0 == matrix_component );
	    }
	}
    }
}

//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( OperatorTools, spectral_radius_test )
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
    matrix_values[0] = -1.0;
    matrix_values[1] = 2.0;
    matrix_values[2] = -1.0;
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

    // Build the Jacobi split.
    Teuchos::RCP<Chimera::LinearOperatorSplit<double,int,int> > a_split =
	Teuchos::rcp( new Chimera::JacobiSplit<double,int,int>( A ) );

    // Check the spectral radius of the iteration matrix.
    double spectral_radius = 
	Chimera::OperatorTools::spectralRadius( a_split->iterationMatrix() );
    TEST_ASSERT( 0.0 <= spectral_radius );
}

//---------------------------------------------------------------------------//
// end tstOperatorTools.cpp
//---------------------------------------------------------------------------//

