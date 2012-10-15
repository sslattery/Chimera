//---------------------------------------------------------------------------//
/*!
 * \file tstBoostRNG.cpp
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

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( BoostRNG, boost_mt11213b_test )
{
    using namespace Chimera;

    boost::mt11213b generator_1;
    Teuchos::RCP<boost::mt11213b> generator_2 = 
	RNGTraits<boost::mt11213b>::create();

    TEST_ASSERT( generator_1.min() == 
		 RNGTraits<boost::mt11213b>::min( *generator_2 ) );

    TEST_ASSERT( generator_1.max() == 
		 RNGTraits<boost::mt11213b>::max( *generator_2 ) );

    RNGTraits<boost::mt11213b>::result_type first_val_1 = generator_1();
    RNGTraits<boost::mt11213b>::result_type first_val_2 = 
	RNGTraits<boost::mt11213b>::generate( *generator_2 );

    for ( int i = 0; i < 10000; ++i )
    {
	TEST_ASSERT( generator_1() == 
		     RNGTraits<boost::mt11213b>::generate( *generator_2 ) );
    }

    RNGTraits<boost::mt11213b>::result_type new_seed = 184839;
    generator_1.seed( new_seed );
    RNGTraits<boost::mt11213b>::setSeed( *generator_2, new_seed );

    RNGTraits<boost::mt11213b>::result_type second_val_1 = generator_1();
    RNGTraits<boost::mt11213b>::result_type second_val_2 = 
	RNGTraits<boost::mt11213b>::generate( *generator_2 );

    TEST_ASSERT( first_val_1 != second_val_1 );
    TEST_ASSERT( first_val_2 != second_val_2 );

    for ( int i = 0; i < 10000; ++i )
    {
	TEST_ASSERT( generator_1() == 
		     RNGTraits<boost::mt11213b>::generate( *generator_2 ) );
    }
}

TEUCHOS_UNIT_TEST( BoostRNG, boost_mt19937_test )
{
    using namespace Chimera;

    boost::mt19937 generator_1;
    Teuchos::RCP<boost::mt19937> generator_2 = 
	RNGTraits<boost::mt19937>::create();

    TEST_ASSERT( generator_1.min() == 
		 RNGTraits<boost::mt19937>::min( *generator_2 ) );

    TEST_ASSERT( generator_1.max() == 
		 RNGTraits<boost::mt19937>::max( *generator_2 ) );

    RNGTraits<boost::mt19937>::result_type first_val_1 = generator_1();
    RNGTraits<boost::mt19937>::result_type first_val_2 = 
	RNGTraits<boost::mt19937>::generate( *generator_2 );

    for ( int i = 0; i < 10000; ++i )
    {
	TEST_ASSERT( generator_1() == 
		     RNGTraits<boost::mt19937>::generate( *generator_2 ) );
    }

    RNGTraits<boost::mt19937>::result_type new_seed = 184839;
    generator_1.seed( new_seed );
    RNGTraits<boost::mt19937>::setSeed( *generator_2, new_seed );

    RNGTraits<boost::mt19937>::result_type second_val_1 = generator_1();
    RNGTraits<boost::mt19937>::result_type second_val_2 = 
	RNGTraits<boost::mt19937>::generate( *generator_2 );

    TEST_ASSERT( first_val_1 != second_val_1 );
    TEST_ASSERT( first_val_2 != second_val_2 );

    for ( int i = 0; i < 10000; ++i )
    {
	TEST_ASSERT( generator_1() == 
		     RNGTraits<boost::mt19937>::generate( *generator_2 ) );
    }
}

//---------------------------------------------------------------------------//
// end tstBoostRNG.cpp
//---------------------------------------------------------------------------//

