//---------------------------------------------------------------------------//
/*!
 * \file tstHistoryBank.cpp
 * \author Stuart R. Slattery
 * \brief History bank tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <Chimera_HistoryBank.hpp>
#include <Chimera_History.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_TypeTraits.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_as.hpp>

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( HistoryBank, history_bank_test )
{
    using namespace Chimera;

    typedef History<double,int> HistoryType;

    int num_histories = 2;

    HistoryBank<HistoryType> bank;

    // Test the bank with an array of histories.
    Teuchos::Array<HistoryType> histories( num_histories );
    for ( int i = 0; i < num_histories; ++i )
    {
	histories[i].setWeight( i );
	histories[i].setGlobalState( i );
    }

    bank.setStack( histories );

    histories.clear();

    TEST_ASSERT( !bank.empty() );
    TEST_ASSERT( Teuchos::as<int>(bank.size()) == num_histories );

    HistoryType history;
    int rindex;
    for ( int i = 0; i < num_histories; ++i )
    {
	rindex = num_histories - i - 1;

	TEST_ASSERT( !bank.empty() );

	TEST_ASSERT( bank.top().weight() == 
		     Teuchos::as<HistoryType::scalar_type>(rindex) );
	TEST_ASSERT( bank.top().globalState() == 
		     Teuchos::as<HistoryType::global_ordinal_type>(rindex) );

	history = bank.pop();
	TEST_ASSERT( history.weight() == 
		     Teuchos::as<HistoryType::scalar_type>(rindex) );
	TEST_ASSERT( history.globalState() == 
		     Teuchos::as<HistoryType::global_ordinal_type>(rindex) );

	TEST_ASSERT( Teuchos::as<int>(bank.size()) == rindex );
    }

    TEST_ASSERT( bank.empty() );
    TEST_ASSERT( bank.size() == 0 );

    // Now test the bank by pushing histories.
    for ( int i = 0; i < num_histories; ++i )
    {
	bank.push( HistoryType( i, i ) );
    }

    TEST_ASSERT( !bank.empty() );
    TEST_ASSERT( Teuchos::as<int>(bank.size()) == num_histories );

    for ( int i = 0; i < num_histories; ++i )
    {
	rindex = num_histories - i - 1;

	TEST_ASSERT( !bank.empty() );

	TEST_ASSERT( bank.top().weight() == 
		     Teuchos::as<HistoryType::scalar_type>(rindex) );
	TEST_ASSERT( bank.top().globalState() == 
		     Teuchos::as<HistoryType::global_ordinal_type>(rindex) );

	history = bank.pop();
	TEST_ASSERT( history.weight() == 
		     Teuchos::as<HistoryType::scalar_type>(rindex) );
	TEST_ASSERT( history.globalState() == 
		     Teuchos::as<HistoryType::global_ordinal_type>(rindex) );

	TEST_ASSERT( Teuchos::as<int>(bank.size()) == rindex );
    }

    TEST_ASSERT( bank.empty() );
    TEST_ASSERT( bank.size() == 0 );

    // Finally test the bank with the stack constructor.
    histories.resize( num_histories );
    for ( int i = 0; i < num_histories; ++i )
    {
	histories[i].setWeight( i );
	histories[i].setGlobalState( i );
    }

    HistoryBank<HistoryType> bank2( histories );

    histories.clear();

    TEST_ASSERT( !bank2.empty() );
    TEST_ASSERT( Teuchos::as<int>(bank2.size()) == num_histories );

    for ( int i = 0; i < num_histories; ++i )
    {
	rindex = num_histories - i - 1;

	TEST_ASSERT( !bank2.empty() );

	TEST_ASSERT( bank2.top().weight() == 
		     Teuchos::as<HistoryType::scalar_type>(rindex) );
	TEST_ASSERT( bank2.top().globalState() == 
		     Teuchos::as<HistoryType::global_ordinal_type>(rindex) );

	history = bank2.pop();
	TEST_ASSERT( history.weight() == 
		     Teuchos::as<HistoryType::scalar_type>(rindex) );
	TEST_ASSERT( history.globalState() == 
		     Teuchos::as<HistoryType::global_ordinal_type>(rindex) );

	TEST_ASSERT( Teuchos::as<int>(bank2.size()) == rindex );
    }

    TEST_ASSERT( bank2.empty() );
    TEST_ASSERT( bank2.size() == 0 );
}

//---------------------------------------------------------------------------//
// end tstHistoryBank.cpp
//---------------------------------------------------------------------------//
