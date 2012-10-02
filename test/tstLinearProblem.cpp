//---------------------------------------------------------------------------//
/*!
 * \file tstLinearProblem.cpp
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

#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( LinearProblem, linear_problem_test )
{
    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();

    // Setup linear system distribution.
    int local_num_rows = 10;
    Teuchos::RCP<const Tpetra::Map<int> > row_map = 
	Tpetra::createLocalMap<int,int>( local_num_rows, comm );

    // Build the identity matrix operator.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int> > A = 
	Tpetra::createCrsMatrix<double,int>( row_map, 1 );
    Teuchos::Array<int> global_columns(1);
    Teuchos::Array<double> matrix_values(1, 1.0);
    int gid = 0;
    for ( int i = 0; i < local_num_rows; ++i )
    {
	gid = i + local_num_rows*comm_rank;
	global_columns[0] = gid;
	A->insertGlobalValues( gid, global_columns(), matrix_values() );
    }
    A->fillComplete();

    // Build the solution vector.
    double X_val = 2.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > X = 
	Tpetra::createVector<double,int>( row_map );
    X->putScalar( X_val );

    // Build the right hand side.
    double B_val = 5.0;
    Teuchos::RCP<Tpetra::Vector<double,int> > B = 
	Tpetra::createVector<double,int>( row_map );
    B->putScalar( B_val );

    // Build the linear problem.
    Teuchos::RCP<Chimera::LinearProblem<double,int> > linear_problem = 
	Teuchos::rcp( new Chimera::LinearProblem<double,int>( A, X, B ) );

    // Check the operator.
    Teuchos::RCP<Tpetra::CrsMatrix<double,int> > A_from_LP = 
	linear_problem->getOperator();
    Tpetra::Vector<double,int> A_diag( row_map );
    A_from_LP->getLocalDiagCopy( A_diag );
    Teuchos::ArrayRCP<const double> diag_view = A_diag.get1dView();
    for ( int i = 0; i < local_num_rows; ++i )
    {
	TEST_ASSERT( diag_view[i] == 1.0 );
    }

    // Check the solution vector.
    Teuchos::RCP<Tpetra::Vector<double,int> > X_from_LP = 
	linear_problem->getLHS();
    Teuchos::ArrayRCP<const double> X_view = X_from_LP->get1dView();
    for ( int i = 0; i < local_num_rows; ++i )
    {
	TEST_ASSERT( X_view[i] == X_val );
    }

    // Check the right-hand side.
    Teuchos::RCP<Tpetra::Vector<double,int> > B_from_LP = 
	linear_problem->getRHS();
    Teuchos::ArrayRCP<const double> B_view = B_from_LP->get1dView();
    for ( int i = 0; i < local_num_rows; ++i )
    {
	TEST_ASSERT( B_view[i] == B_val );
    }

    // Check that the residual initialized to zero.
    Teuchos::RCP<Tpetra::Vector<double,int> > R_from_LP = 
	linear_problem->getResidual();
    Teuchos::ArrayRCP<const double> R_view = R_from_LP->get1dView();
    for ( int i = 0; i < local_num_rows; ++i )
    {
	TEST_ASSERT( R_view[i] == 0.0 );
    }
    
    // Compute the residual.
    linear_problem->computeResidual();

    // Check the residual.
    double R_val = B_val - X_val;
    Teuchos::RCP<Tpetra::Vector<double,int> > new_R_from_LP = 
	linear_problem->getResidual();
    Teuchos::ArrayRCP<const double> new_R_view = new_R_from_LP->get1dView();
    for ( int i = 0; i < local_num_rows; ++i )
    {
	TEST_ASSERT( new_R_view[i] == R_val );
    }
}

//---------------------------------------------------------------------------//
// end tstLinearProblem.cpp
//---------------------------------------------------------------------------//

