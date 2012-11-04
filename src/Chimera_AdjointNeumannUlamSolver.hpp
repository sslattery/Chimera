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
 * \file Chimera_AdjointNeumannUlamSolver.hpp
 * \author Stuart R. Slattery
 * \brief Adjoint Neumann-Ulam solver declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_ADJOINTNEUMANNULAMSOLVER_HPP
#define Chimera_ADJOINTNEUMANNULAMSOLVER_HPP

#include "Chimera_NeumannUlamSolver.hpp"
#include "Chimera_LinearProblem.hpp"
#include "Chimera_LinearOperatorSplit.hpp"
#include "Chimera_History.hpp"
#include "Chimera_HistoryBank.hpp"
#include "Chimera_HistoryBuffer.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_CrsMatrix.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class NeumannUlamSolver
 * \brief Interface definition for Neumann-Ulam Monte Carlo solvers.
 */
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO, class RNG>
class AdjointNeumannUlamSolver : public NeumannUlamSolver<Scalar,LO,GO,RNG>
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                    scalar_type;
    typedef LO                                        local_ordinal_type;
    typedef GO                                        global_ordinal_type;
    typedef RNG                                       rng_type;
    typedef NeumannUlamSolver<Scalar,LO,GO,RNG>       Base;
    typedef typename Base::RCP_LinearProblem          RCP_LinearProblem;
    typedef typename Base::RCP_LinearOperatorSplit    RCP_LinearOperatorSplit;
    typedef typename Base::RCP_RNG                    RCP_RNG;
    typedef Teuchos::RCP<Teuchos::ParameterList>      RCP_ParameterList;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO>           TpetraCrsMatrix;
    typedef Teuchos::RCP<TpetraCrsMatrix>             RCP_TpetraCrsMatrix;
    typedef History<Scalar,GO>                        HistoryType;
    //@}

    //! Constructor.
    AdjointNeumannUlamSolver( 
	const RCP_LinearProblem& linear_problem,
	const RCP_LinearOperatorSplit& linear_operator_split,
	const RCP_RNG& rng,
	const RCP_ParameterList& plist );

    //! Destructor.
    ~AdjointNeumannUlamSolver();

    //! Execute a stage of histories with random walks.
    void walk();

  private:

    // Build the probability matrix.
    void buildProbabilityMatrix();

    // Build the ghosted iteration matrix.
    void buildGhostIterationMatrix();

    // Sample the source to build a starting history bank.
    HistoryBank<HistoryType> sampleSource();

    // Check for completion of all random walks.
    bool allBanksEmpty( const HistoryBank<HistoryType>& bank );

    // Check for empty buffers on all processes.
    bool allBuffersEmpty( const HistoryBuffer<HistoryType>& buffer );

  private:

    // Relative weight cutoff.
    Scalar d_relative_weight_cutoff;

    // Probability matrix.
    RCP_TpetraCrsMatrix d_probability_matrix;

    // Ghosted iteration matrix.
    RCP_TpetraCrsMatrix d_ghost_iteration_matrix;
};

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_AdjointNeumannUlamSolver_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_ADJOINTNEUMANNULAMSOLVER_HPP

//---------------------------------------------------------------------------//
// end Chimera_AdjointNeumannUlamSolver.hpp
//---------------------------------------------------------------------------//
