//----------------------------------*-C++-*----------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
/*!
 * \file   mesh/test/tstEpetra.cpp
 * \author Stuart Slattery
 * \brief  Epetra class unit tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

TEUCHOS_UNIT_TEST( Epetra_Vector, epetra_vector_test)
{
    int error = 0;

    Teuchos::ArrayRCP<double> data(10, 1.3);
    Epetra_SerialComm comm;
    Epetra_Map map( (int) data.size(), 0, comm );
    Epetra_Vector epetra_vector( View, map, data.get() );

    Teuchos::ArrayRCP<double> update_array(10, 1.0);
    Epetra_Vector update_vector( View, map, update_array.get() );

    Teuchos::ArrayRCP<double>::const_iterator const_it;
    for ( int i = 0; i < (int) data.size(); ++i )
    {
	TEST_ASSERT( epetra_vector[i] == 1.3 );
    }

    std::vector<double> values( 4, 2.5 );
    std::vector<int> indices( 4 );
    indices[0] = 3;
    indices[1] = 4;
    indices[2] = 5;
    indices[3] = 6;
    error = epetra_vector.ReplaceGlobalValues( 4, 0, &values[0], &indices[0] );
    TEST_ASSERT( error == 0 );
    for ( int i = 3; i < 7; ++i )
    {
	TEST_ASSERT( epetra_vector[i] == 2.5 );
	TEST_ASSERT( data[i] == 2.5 );
    }

    error = epetra_vector.SumIntoGlobalValues( 4, 0, &values[0], &indices[0] );
    TEST_ASSERT( error == 0 );
    for ( int i = 3; i < 7; ++i )
    {
	TEST_ASSERT( epetra_vector[i] == 5.0 );
	TEST_ASSERT( data[i] == 5.0 );
    }

    error = epetra_vector.Scale( 2.0 );
    TEST_ASSERT( error == 0 );
    for ( int i = 0; i < 3; ++i )
    {
	TEST_ASSERT( epetra_vector[i] == 2.6 );
	TEST_ASSERT( data[i] == 2.6 );
	TEST_ASSERT( epetra_vector[i+7] == 2.6 );
	TEST_ASSERT( data[i+7] == 2.6 );

    }
    for ( int i = 3; i < 7; ++i )
    {
	TEST_ASSERT( epetra_vector[i] == 10.0 );
	TEST_ASSERT( data[i] == 10.0 );
    }

    error = epetra_vector.Update( -1.0, update_vector, 1.0 );
    TEST_ASSERT( error == 0 );
    for ( int i = 3; i < 7; ++i )
    {
	TEST_ASSERT( epetra_vector[i] == 9.0 );
	TEST_ASSERT( data[i] == 9.0 );
    }

    double norm_2 = 0.0;
    for ( int i = 0; i < (int) data.size(); ++i )
    {
	norm_2 += epetra_vector[i]*epetra_vector[i];
    }
    norm_2 = pow( norm_2, 0.5 );
    double l2_norm = 0;
    error = epetra_vector.Norm2( &l2_norm );
    TEST_ASSERT( error == 0 );
    TEST_ASSERT( l2_norm == norm_2 );


    double inf_norm = 0;
    error = epetra_vector.NormInf( &inf_norm );
    TEST_ASSERT( error == 0 );
    TEST_ASSERT( inf_norm == 9.0 );
}

TEUCHOS_UNIT_TEST( Epetra_CrsMatrix, epetra_crsmatrix_test)
{
    int error = 0;

    int vector_size = 10;

    Epetra_SerialComm comm;
    Epetra_Map map( vector_size, 0, comm );

    Teuchos::ArrayRCP<double> u_array( vector_size, 2.0 );
    Epetra_Vector u( View, map, u_array.get() );

    std::vector<int> entries_per_row( vector_size, 1 );
    Epetra_CrsMatrix A( Copy, map, &entries_per_row[0] );
    double diag_value = 2.0;
    for ( int i = 0; i < vector_size; ++i )
    {
	error = A.InsertGlobalValues( i, 1, &diag_value, &i );
	TEST_ASSERT( error == 0 );
    }
    error = A.FillComplete();
    TEST_ASSERT( error == 0 );

    Teuchos::ArrayRCP<double> v_array( vector_size );
    Epetra_Vector v( View, map, v_array.get() );
    error = A.Apply( u, v );
    TEST_ASSERT( error == 0 );
    for ( int i = 0; i < vector_size; ++i )
    {
	TEST_ASSERT( v_array[i] == 4.0 );
    }

}

//---------------------------------------------------------------------------//
//                        end of tstEpetra.cpp
//---------------------------------------------------------------------------//
