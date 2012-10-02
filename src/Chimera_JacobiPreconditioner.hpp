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
// \file Chimera_JacobiPreconditioner.hpp
// \author Stuart R. Slattery
// \brief Jacobi preconditioner declaration.
//---------------------------------------------------------------------------//

#ifndef Chimera_JACOBIPRECONDITIONER_HPP
#define Chimera_JACOBIPRECONDITIONER_HPP

#include <Teuchos_RCP.hpp>

#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

namespace Chimera
{
//---------------------------------------------------------------------------//
// 
//---------------------------------------------------------------------------//
class JacobiPreconditioner
{
  private:

    // Linear problem.
    Teuchos::RCP<Epetra_LinearProblem> d_linear_problem;

    // Preconditioned operator.
    Teuchos::RCP<Epetra_CrsMatrix> d_M_inv_A;

    // Preconditioned rhs.
    Teuchos::RCP<Epetra_Vector> d_M_inv_b;

    // Operator diagonal elements.
    Teuchos::RCP<Epetra_Vector> d_diagonal;

  public:

    // Constructor.
    JacobiPreconditioner( Teuchos::RCP<Epetra_LinearProblem> &linear_problem );

    // Destructor.
    ~JacobiPreconditioner();

    // Precondition both the operator and the right hand side.
    void precondition();

    // Do preconditioning on the operator.
    void preconditionOperator();

    // Do preconditioning on the right hand side.
    void preconditionRHS();

    // Get the preconditioned operator.
    Teuchos::RCP<Epetra_CrsMatrix> getOperator() const
    { return d_M_inv_A; }

    // Get the preconditioned rhs.
    Teuchos::RCP<Epetra_Vector> getRHS() const
    { return d_M_inv_b; }
};

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_JACOBIPRECONDITIONER_HPP

//---------------------------------------------------------------------------//
// end Chimera_JacobiPreconditioner.hpp
//---------------------------------------------------------------------------//


