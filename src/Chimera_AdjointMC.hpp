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
// \file Chimera_AdjointMC.hpp
// \author Stuart Slattery
// \brief Adjoint Monte Carlo solver declaration.
//---------------------------------------------------------------------------//

#ifndef Chimera_ADJOINTMC_HPP
#define Chimera_ADJOINTMC_HPP

#include "Chimera_BoostRNG.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

namespace Chimera
{

class AdjointMC
{
  public:

    // Constructor.
    AdjointMC( Teuchos::RCP<Epetra_LinearProblem> &linear_problem,
	       Teuchos::RCP<Teuchos::ParameterList> &plist );

    // Destructor.
    ~AdjointMC();

    // Solve.
    void walk();

    // Return the iteration matrix.
    const Teuchos::RCP<Epetra_CrsMatrix>& getH() const
    { return d_H; }

  private:

    // Build the iteration matrix.
    Teuchos::RCP<Epetra_CrsMatrix> buildH();

    // Build the adjoint probability matrix.
    Epetra_CrsMatrix buildQ();

    // Build the cumulative distribution function.
    Epetra_CrsMatrix buildC();

  private:

    // Linear problem.
    Teuchos::RCP<Epetra_LinearProblem> d_linear_problem;

    // Parameter list.
    Teuchos::RCP<Teuchos::ParameterList> d_plist;

    // Random number generator.
    Teuchos::RCP<boost::mt11213b> d_rng;

    // Iteration matrix.
    Teuchos::RCP<Epetra_CrsMatrix> d_H;

    // Adjoint probability matrix.
    Epetra_CrsMatrix d_Q;

    // Cumulative distribution function.
    Epetra_CrsMatrix d_C;
};

} // end namespace Chimera

#endif // end Chimera_ADJOINTMC_HPP

//---------------------------------------------------------------------------//
// end Chimera_AdjointMC.hpp
//---------------------------------------------------------------------------//

