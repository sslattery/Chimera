//---------------------------------------------------------------------------//
/*
  Copyright (c) 2012, Stuart R. Slattery
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

  *: Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  *: Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  *: Neither the name of the University of Wisconsin - Madison nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//---------------------------------------------------------------------------//
/*!
 * \file Chimera_LinearSolverFactory_def.hpp
 * \author Stuart Slattery
 * \brief Linear solver factory definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_LINEARSOLVERFACTORY_DEF_HPP
#define Chimera_LINEARSOLVERFACTORY_DEF_HPP

#include <string>

#include "Chimera_Assertion.hpp"
#include "Chimera_StationarySolver.hpp"

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Creation method.
 */
template<class Scalar, class LO, class GO>
Teuchos::RCP<LinearSolver<Scalar,LO,GO> >
LinearSolverFactory::create( 
    const Teuchos::RCP<Teuchos::ParameterList>& plist,
    const Teuchos::RCP<LinearProblem<Scalar,LO,GO> >& linear_problem )
{
    testPrecondition( !plist.is_null() );
    testPrecondition( !linear_problem.is_null() );

    Teuchos::RCP<LinearSolver<Scalar,LO,GO> > linear_solver;

    if( plist->get<std::string>("SOLVER TYPE") == "STATIONARY" )
    {
	linear_solver = 
	    Teuchos::rcp( new StationarySolver<Scalar,LO,GO>( linear_problem,
							      plist ) );
    }

    testPostcondition( !linear_solver.is_null() );

    return linear_solver;
}

//---------------------------------------------------------------------------//

} // end namepsace Chimera

#endif // end Chimera_LINEARSOLVERFACTORY_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_LinearSolverFactory_def.hpp
//---------------------------------------------------------------------------//
