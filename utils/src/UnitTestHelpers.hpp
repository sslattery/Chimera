//---------------------------------------------------------------------------//
/*! 
 * \file UnitTestHelpers.hpp
 * \author Stuart R. Slattery
 * \bref Standard unit test #helpers
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_UNITESTHELPERS_HPP
#define CHIMERA_UNITESTHELPERS_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <ostream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>

namespace Chimera
{

namespace UnitTestHelpers
{

//@{
//! typedefs.
typedef Teuchos::RCP<const Teuchos::Comm<int> >       RCP_Comm;
//@}

/*!
 * \brief Get the default communicator.
 */
RCP_Comm getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<int>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<int>() );
#endif
}

} // end namespace UnitTestHelpers

} // end namespace Chimera

#endif // end CHIMERA_UNITESTHELPERS_HPP

//---------------------------------------------------------------------------//
// end UnitTestHelpers.hpp
//---------------------------------------------------------------------------//

