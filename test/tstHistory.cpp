//---------------------------------------------------------------------------//
/*!
 * \file tstHistory.cpp
 * \author Stuart R. Slattery
 * \brief History tests.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <cassert>

#include <Chimera_History.hpp>
#include <Chimera_RNGTraits.hpp>
#include <Chimera_BoostRNG.hpp>

#include <boost/random/mersenne_twister.hpp>

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

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( History, history_test )
{
    using namespace Chimera;

    Teuchos::RCP<boost::mt19937> rng = RNGTraits<boost::mt19937>::create();

    int num_rand = 10;
    Teuchos::Array<double> randoms(num_rand);
    for ( int i = 0; i < num_rand; ++i )
    {
	randoms[i] = 
	    Teuchos::as<double>(RNGTraits<boost::mt19937>::generate(*rng)) /
	    Teuchos::as<double>(RNGTraits<boost::mt19937>::max(*rng));
    }

    History<double,int,int> history;
    TEST_ASSERT( history.active() );

    double total = 2.0;    
    history.setWeight( total );
    TEST_ASSERT( history.weight() == total );

    for ( int i = 0; i < num_rand; ++i )
    {
	history.addWeight( randoms[i] );
	total += randoms[i];
	TEST_ASSERT( history.weight() == total );
    }

    total = 2.0;
    history.setWeight( total );
    for ( int i = 0; i < num_rand; ++i )
    {
	history.multiplyWeight( randoms[i] );
	total *= randoms[i];
	TEST_ASSERT( history.weight() == total );
    }

    int local_state = 2929;
    history.setLocalState( local_state );
    TEST_ASSERT( history.localState() == local_state );

    int global_state = 5943;
    history.setGlobalState( global_state );
    TEST_ASSERT( history.globalState() == global_state );

    history.terminate();
    TEST_ASSERT( !history.active() );
}

//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( History, serialization_test )
{
    using namespace Chimera;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();
    int comm_rank = comm->getRank();
    
    History<double,int,int> history;

    double weight = 4.393939;
    int local_state = 4941;
    int global_state = 28390;

    if ( comm_rank == 0 )
    {
	history.setWeight( weight );
	history.setLocalState( local_state );
	history.setGlobalState( global_state );
    }
    comm->barrier();
    
    Teuchos::broadcast( *comm, 0, 
			Teuchos::Ptr<History<double,int,int> >( &history ) );
 
    TEST_ASSERT( history.weight() == weight );
    TEST_ASSERT( history.localState() == local_state );
    TEST_ASSERT( history.globalState() == global_state );
}

//---------------------------------------------------------------------------//
// end tstHistory.cpp
//---------------------------------------------------------------------------//
