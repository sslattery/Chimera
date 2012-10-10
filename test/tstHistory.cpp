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

//---------------------------------------------------------------------------//
// Tests
//---------------------------------------------------------------------------//
TEUCHOS_UNIT_TEST( History, history_test )
{
    using namespace Chimera;

}

//---------------------------------------------------------------------------//
// end tstHistory.cpp
//---------------------------------------------------------------------------//
