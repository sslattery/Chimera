//---------------------------------------------------------------------------//
/*!
 * \file   Exception.hpp
 * \author Stuart R. Slattery
 * \brief  Exceptions for error handling.
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_EXCEPTION_HPP
#define CHIMERA_EXCEPTION_HPP

#include <stdexcept>
#include <string>

namespace Chimera
{

//---------------------------------------------------------------------------//
/*!
 * \brief Exception class to be thrown when function preconditions are not
 * met.
 */
class PreconditionException : public std::runtime_error
{
  public:
    PreconditionException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

/*!
 * \brief Exception class to be thrown when function postconditions are not
 * met. 
 */
class PostconditionException : public std::runtime_error
{
  public:
    PostconditionException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

/*!
 * \brief Exception class to be thrown when a function alters an invariant.
 */
class InvariantException : public std::runtime_error
{
  public:
    InvariantException( const std::string &msg )
	: std::runtime_error( msg )
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// Test for a precondition exception.
void testPrecondition( bool throw_if_false, const std::string &msg );

// Test for a postcondition exception.
void testPostcondition( bool throw_if_false, const std::string &msg );

// Test for a Invariant exception.
void testInvariant( bool throw_if_false, const std::string &msg );

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end CHIMERA_EXCEPTION_HPP
