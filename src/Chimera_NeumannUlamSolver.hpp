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
 * \file Chimera_NeumannUlamSolver.hpp
 * \author Stuart R. Slattery
 * \brief Neumann-Ulam solver interface definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_NEUMANNULAMSOLVER_HPP
#define Chimera_NEUMANNULAMSOLVER_HPP

#include <Chimera_LinearProblem.hpp>
#include <Chimera_LinearOperatorSplit.hpp>

#include <Teuchos_RCP.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class NeumannUlamSolver
 * \brief Interface definition for Neumann-Ulam Monte Carlo solvers.
 */
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO, class RNG>
class NeumannUlamSolver
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                    scalar_type;
    typedef LO                                        local_ordinal_type;
    typedef GO                                        global_ordinal_type;
    typedef RNG                                       rng_type;
    typedef Teuchos::RCP<RNG>                         RCP_RNG;
    typedef LinearProblem<Scalar,LO,GO>               LinearProblemType;
    typedef Teuchos::RCP<LinearProblemType>           RCP_LinearProblem;
    typedef LinearOperatorSplit<Scalar,LO,GO>         LinearOperatorSplitType;
    typedef Teuchos::RCP<LinearOperatorSplitType>     RCP_LinearOperatorSplit;
    //@}

    //! Constructor.
    NeumannUlamSolver()
    { /* ... */ }

    //! Destructor.
    virtual ~NeumannUlamSolver()
    { /* ... */ }

    //! Execute a stage of histories with random walks.
    virtual void walk() = 0;

    //! Get the random number generator.
    RCP_RNG rng() const
    { return b_rng; }

    //! Get the linear problem.
    RCP_LinearProblem linearProblem() const
    { return b_linear_problem; }

    //! Get the linear operator split.
    RCP_LinearOperatorSplit linearOperatorSplit() const
    { return b_linear_operator_split; }

    //! Get the weight cutoff
    Scalar weightCutoff() const
    { return b_weight_cutoff; }

    //! Get the number of histories per stage.
    GO historiesPerStage() const
    { return b_histories_per_stage; }

  protected:

    // The linear problem.
    RCP_LinearProblem b_linear_problem;

    // Linear operator split.
    RCP_LinearOperatorSplit b_linear_operator_split;

    // Random number generator.
    RCP_RNG b_rng;

    // Cutoff weight.
    Scalar b_weight_cutoff;

    // Number of histories per stage.
    int b_histories_per_stage;
};

} // end namespace Chimera

#endif // end Chimera_NEUMANNULAMSOLVER_HPP

//---------------------------------------------------------------------------//
// end Chimera_NeumannUlamSolver.hpp
//---------------------------------------------------------------------------//

