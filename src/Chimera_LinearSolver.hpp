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
 * \file Chimera_LinearSolver.hpp
 * \author Stuart R. Slattery
 * \brief LinearSolver interface definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_LINEARSOLVER_HPP
#define Chimera_LINEARSOLVER_HPP

#include <Chimera_LinearProblem.hpp>

#include <Teuchos_RCP.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class LinearSolver
 * \brief Interface definition for sparse iterative linear solvers.
 */
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO>
class LinearSolver
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                  scalar_type;
    typedef LO                                      local_ordinal_type;
    typedef GO                                      global_ordinal_type;
    typedef LinearProblem<Scalar,LO,GO>             LinearProblemType;
    typedef Teuchos::RCP<LinearProblemType>         RCP_LinearProblem;
    //@}

    //! Constructor.
    LinearSolver()
    { /* ... */ }

    //! Destructor.
    virtual ~LinearSolver()
    { /* ... */ }

    //! Set the convergence tolerance
    void setTolerance( Scalar tolerance )
    { b_tolerance = tolerance; }

    //! Iterate until convergence.
    virtual void iterate() = 0;

    //! Get the linear problem.
    RCP_LinearProblem linearProblem() const
    { return b_linear_problem; }

    //! Get the convergence tolerance.
    Scalar tolerance() const
    { return b_tolerance; }

    //! Get the number of iterations needed to converge.
    int numIters() const
    { return b_num_iters; }

  protected:

    // The linear problem.
    RCP_LinearProblem b_linear_problem;

    // Convergence tolerance.
    Scalar b_tolerance;

    // Number of iterations to converge last solution.
    int b_num_iters;
};

} // end namespace Chimera

#endif // end Chimera_LINEARSOLVER_HPP

//---------------------------------------------------------------------------//
// end Chimera_LinearSolver.hpp
//---------------------------------------------------------------------------//

