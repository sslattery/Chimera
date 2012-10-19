//---------------------------------------------------------------------------//
/*!
 * \file tstHistoryBuffer.cpp
 * \author Stuart R. Slattery
 * \brief History buffer tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <Chimera_HistoryBuffer.hpp>
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
TEUCHOS_UNIT_TEST( HistoryBuffer, history_buffer_test )
{
    using namespace Chimera;

    typedef History<double,int> HistoryType;

    // Setup parallel distribution.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    int comm_size = comm->getSize();

    // Setup map.
    int local_num_histories = 10;
    int global_num_histories = local_num_histories*comm_size;
    Teuchos::RCP<const Tpetra::Map<int> > map = 
	Tpetra::createUniformContigMap<int,int>( global_num_histories, comm );

    // Build a history buffer with inverse rank destination states.
    HistoryBuffer<HistoryType> history_buffer();
    int inverse_rank = comm_size - comm_rank - 1;
    double weight = comm_rank*1.0;
    int global_state = inverse_rank*local_num_histories;
    for ( int i = 0; i < local_num_histories; ++i )
    {
	history_buffer.pushBack( 
	    HistoryType( weight, global_state ) );
    }

    // Communicate the histories to their destinations. 
    HistoryBank<HistoryType> bank = history_buffer.communicate( map );

    // Check the bank created.
    TEST_ASSERT( !bank.empty() );
    TEST_ASSERT( Teuchos::as<int>(bank.size()) == local_num_histories );
    HistoryType history;
    for ( int i = 0; i < local_num_histories; ++i )
    {
	TEST_ASSERT( !bank.empty() );

	history = bank.pop();
	TEST_ASSERT( history.weight() == 
		     Teuchos::as<HistoryType::scalar_type>(inverse_rank) );
	TEST_ASSERT( history.globalState() == 
		     Teuchos::as<HistoryType::global_ordinal_type>(
			 comm_rank*local_num_histories) );
    }

    TEST_ASSERT( bank.empty() );
}

//---------------------------------------------------------------------------//
// end tstHistoryBuffer.cpp
//---------------------------------------------------------------------------//
