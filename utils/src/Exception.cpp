//---------------------------------------------------------------------------//
/*!
 * \file   Exception.cpp
 * \author Stuart R. Slattery
 * \brief  Exceptions for error handling.
 */
//---------------------------------------------------------------------------//

#include <iostream>

#include "Exception.hpp"

#include <Teuchos_TestForException.hpp>

namespace Chimera
{
/*!
 * \brief Test for a precondition exception.
 */
void testPrecondition( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				PreconditionException,
				msg << std::endl );
}

/*!
 * \brief Test for a postcondition exception.
 */
void testPostcondition( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				PostconditionException,
				msg << std::endl );
}

/*!
 * \brief Test for a Invariant exception.
 */
void testInvariant( bool throw_if_false, const std::string &msg )
{
    TEUCHOS_TEST_FOR_EXCEPTION( !throw_if_false, 
				InvariantException,
				msg << std::endl );
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end Exception.cpp
//---------------------------------------------------------------------------//
