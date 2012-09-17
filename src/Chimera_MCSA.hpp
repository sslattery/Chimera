//---------------------------------------------------------------------------//
// \file Chimera_MCSA.hpp
// \author Stuart R. Slattery
// \brief Monte Carlo Synthetic Acceleration solver declaration.
//---------------------------------------------------------------------------//

#ifndef Chimera_MCSA_HPP
#define Chimera_MCSA_HPP

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

namespace Chimera
{
namespace Solvers
{

class MCSA
{
  private:

    // Linear problem.
    Teuchos::RCP<Epetra_LinearProblem> d_linear_problem;

    // Iteration count.
    int d_num_iters;

  public:

    // Constructor.
    MCSA( Teuchos::RCP<Epetra_LinearProblem> &linear_problem );

    // Destructor.
    ~MCSA();

    // Solve.
    void iterate( const int max_iters, const double tolerance,
		  const int num_histories, const double weight_cutoff );

    // Get the iteration count from the last solve.
    int getNumIters() const
    { return d_num_iters; }
};

//---------------------------------------------------------------------------//

} // end namespace Solvers
} // end namespace Chimera

#endif // end Chimera_MCSA_HPP

//---------------------------------------------------------------------------//
// end Chimera_MCSA.hpp
//---------------------------------------------------------------------------//

