//---------------------------------------------------------------------------//
/*!
 * \file   Chimera_Assertion.hpp
 * \author Stuart Slattery
 * \brief  Assertions and Design-by-Contract for error handling.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_ASSERTION_HPP
#define Chimera_ASSERTION_HPP

#include <stdexcept>
#include <string>

#include "Chimera_config.hpp"

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Base class for Chimera assertions. This structure is heavily based on
 * that in Nemesis developed by Tom Evans. We derive from std::logic_error
 * here as the DBC checks that utilize this class are meant to find errors
 * that can be prevented before runtime.
 */
//---------------------------------------------------------------------------//
class Assertion : public std::logic_error
{
  public:

    /*! 
     * \brief Default constructor.
     *
     * \param msg Error message.
     */
    Assertion( const std::string& msg )
	: std::logic_error( msg )
    { /* ... */ }

    /*! 
     * \brief DBC constructor.
     *
     * \param cond A string containing the assertion condition that failed.
     *
     * \param field A string containing the file name in which the assertion
     * failed. 
     *
     * \param line The line number at which the assertion failed.
     */
    Assertion( const std::string& cond, const std::string& file, 
	       const int line )
	: std::logic_error( generate_output( cond, file, line ) )
    { /* ... */ }

    //! Destructor.
    virtual ~Assertion() throw()
    { /* ... */ }

  private:

    // Build an assertion output from advanced constructor arguments.
    std::string generate_output( const std::string& cond, 
				 const std::string& file, 
				 const int line ) const;
};

//---------------------------------------------------------------------------//
// Throw functions.
//---------------------------------------------------------------------------//
// Throw a Chimera::Assertion.
void throwAssertion( const std::string& cond, const std::string& file,
		     const int line );

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Design-by-Contract macros.
//---------------------------------------------------------------------------//
/*!
  \page Chimera Design-by-Contract.
 
  Design-by-Contract (DBC) functionality is provided to verify function
  preconditions, postconditions, and invariants. These checks are separated
  from the debug build and can be activated for both release and debug
  builds. They can be activated by setting the following in a CMake
  configure:
 
  -D Chimera_ENABLE_DBC:BOOL=ON
 
  By default, DBC is deactivated. Although they will require additional
  computational overhead, these checks provide a mechanism for veryifing
  library input arguments. Note that the bounds-checking functionality used
  within the Chimera is only provided by a debug build.
 
  In addition, rememberValue is provided to store values used only for DBC
  checks and no other place in executed code.

  Separate from the DBC build, testAssertion can be used at any time verify a
  conditional. This should be used instead of the standard cassert.
 */

#if HAVE_CHIMERA_DBC

#define testPrecondition(c) \
    if (!(c)) Chimera::throwAssertion( #c, __FILE__, __LINE__ )
#define testPostcondition(c) \
    if (!(c)) Chimera::throwAssertion( #c, __FILE__, __LINE__ )
#define testInvariant(c) \
    if (!(c)) Chimera::throwAssertion( #c, __FILE__, __LINE__ )
#define rememberValue(c) c

#else

#define testPrecondition(c)
#define testPostcondition(c)
#define testInvariant(c)
#define rememberValue(c)

#endif


#define testAssertion(c) \
    if (!(c)) Chimera::throwAssertion( #c, __FILE__, __LINE__ )

//---------------------------------------------------------------------------//

#endif // end Chimera_ASSERTION_HPP

//---------------------------------------------------------------------------//
// end Chimera_Assertion.hpp
//---------------------------------------------------------------------------//

