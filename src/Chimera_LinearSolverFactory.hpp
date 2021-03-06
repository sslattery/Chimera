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
 * \file Chimera_LinearSolverFactory.hpp
 * \author Stuart Slattery
 * \brief Linear solver factory declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_LINEARSOLVERFACTORY_HPP
#define Chimera_LINEARSOLVERFACTORY_HPP

#include "Chimera_LinearSolver.hpp"
#include "Chimera_LinearProblem.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class LinerarSolverFactory
 * \brief Factory for linear operator split implementations.
 */
//---------------------------------------------------------------------------//
class LinearSolverFactory
{
  public:

    //! Constructor.
    LinearSolverFactory()
    { /* ... */ }

    //! Destructor.
    ~LinearSolverFactory()
    { /* ... */ }

    // Factory method.
    template<class Scalar, class LO, class GO>
    static Teuchos::RCP<LinearSolver<Scalar,LO,GO> >
    create( const Teuchos::RCP<Teuchos::ParameterList>& plist,
	    const Teuchos::RCP<LinearProblem<Scalar,LO,GO> >& linear_problem );
};

} // end namepsace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_LinearSolverFactory_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_LINEARSOLVERFACTORY_HPP

//---------------------------------------------------------------------------//
// end Chimera_LinearSolverFactory.hpp
//---------------------------------------------------------------------------//


