//---------------------------------------------------------------------------//
// \file JacobiPreconditioner.hpp
// \author Stuart R. Slattery
// \brief Jacobi preconditioner declaration.
//---------------------------------------------------------------------------//

#ifndef HMCSA_JACOBIPRECONDITIONER_HPP
#define HMCSA_JACOBIPRECONDITIONER_HPP

#include <Teuchos_RCP.hpp>

#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

namespace HMCSA
{

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

} // End namespace HMCSA

#endif // end HMCSA_JACOBIPRECONDITIONER_HPP

//---------------------------------------------------------------------------//
// end JacobiPreconditioner.hpp
//---------------------------------------------------------------------------//


