//---------------------------------------------------------------------------//
// \file AdjointMC.hpp
// \author Stuart Slattery
// \brief Adjoint Monte Carlo solver declaration.
//---------------------------------------------------------------------------//

#ifndef HMCSA_ADJOINTMC_HPP
#define HMCSA_ADJOINTMC_HPP

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

namespace HMCSA
{

class AdjointMC
{
  private:

    // Linear problem.
    Teuchos::RCP<Epetra_LinearProblem> d_linear_problem;

    // Iteration matrix.
    Epetra_CrsMatrix d_H;

    // Adjoint probability matrix.
    Epetra_CrsMatrix d_Q;

    // Cumulative distribution function.
    Epetra_CrsMatrix d_C;

  public:

    // Constructor.
    AdjointMC( Teuchos::RCP<Epetra_LinearProblem> &linear_problem );

    // Destructor.
    ~AdjointMC();

    // Solve.
    void walk( const int num_histories, const double weight_cutoff );

    // Return the iteration matrix.
    const Epetra_CrsMatrix& getH() const
    { return d_H; }

  private:

    // Build the iteration matrix.
    Epetra_CrsMatrix buildH();

    // Build the adjoint probability matrix.
    Epetra_CrsMatrix buildQ();

    // Build the cumulative distribution function.
    Epetra_CrsMatrix buildC();
};

} // end namespace HMCSA

#endif // end HMCSA_ADJOINTMC_HPP

//---------------------------------------------------------------------------//
// end AdjointMC.hpp
//---------------------------------------------------------------------------//

