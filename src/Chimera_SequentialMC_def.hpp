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
 * \file Chimera_SequentialMC_def.hpp
 * \author Stuart R. Slattery
 * \brief SequentialMC definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_SEQUENTIALMC_DEF_HPP
#define Chimera_SEQUENTIALMC_DEF_HPP

#include "Chimera_Assertion.hpp"
#include "Chimera_RNGTraits.hpp"
#include "Chimera_LinearOperatorSplitFactory.hpp"
#include "Chimera_NeumannUlamSolverFactory.hpp"

#include <Teuchos_ScalarTraits.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO, class RNG>
SequentialMC<Scalar,LO,GO,RNG>::SequentialMC( 
    const RCP_LinearProblem& linear_problem,
    const RCP_ParameterList& plist )
    : d_rng( RNGTraits<RNG>::create() )
{
    // Check preconditions.
    testPrecondition( !linear_problem.is_null() );
    testPrecondition( !plist.is_null() );

    // Set the base data.
    this->b_linear_problem = linear_problem;
    this->b_linear_operator_split = LinearOperatorSplitFactory::create( 
	plist, this->b_linear_problem->getOperator() );
    this->b_tolerance = plist->get<Scalar>("TOLERANCE");
    this->b_max_num_iters = plist->get<int>("MAX ITERS");
    this->b_num_iters = 0;
    this->b_is_converged = false;

    // Generate the residual linear problem.
    Teuchos::RCP<const Tpetra::Map<LO,GO> > row_map =
	this->b_linear_problem->getOperator()->getRowMap();

    Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> > delta_X =
	Tpetra::createVector<Scalar,LO,GO>( row_map );

    RCP_LinearProblem residual_problem = 
	Teuchos::rcp( new LinearProblem<Scalar,LO,GO>( 
			  this->b_linear_problem->getOperator(),
			  delta_X,
			  this->b_linear_problem->getResidual() ) );

    // Build the Neumann-Ulam solver.
    d_nu_solver = NeumannUlamSolverFactory::create( 
	plist, residual_problem, this->b_linear_operator_split, d_rng );

    // Check postconditions.
    testPostcondition( !this->b_linear_operator_split.is_null() );
    testPostcondition( !d_rng.is_null() );
    testPostcondition( !d_nu_solver.is_null() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Scalar, class LO, class GO, class RNG>
SequentialMC<Scalar,LO,GO,RNG>::~SequentialMC()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Iterate until convergence.
 */
template<class Scalar, class LO, class GO, class RNG>
void SequentialMC<Scalar,LO,GO,RNG>::iterate()
{
    // Generate the convergence criteria.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType b_norm = 
	this->b_linear_problem->getRHS()->normInf();
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
	convergence_criteria = this->b_tolerance * b_norm;

    // Setup for iteration.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
	residual_norm = 1.0;
    this->b_num_iters = 0;
    this->b_is_converged = false;

    // Iterate.
    while ( residual_norm > convergence_criteria &&
	    this->b_num_iters < this->b_max_num_iters )
    {
	// Compute the new residual.
	Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> > res_vec = Teuchos::rcp(
	    new Tpetra::Vector<Scalar,LO,GO>( *this->b_linear_problem->getResidual() ) );
	this->b_linear_problem->computeResidual();
	this->b_linear_operator_split->applyInvM(
	    this->b_linear_problem->getResidual(),
	    this->b_linear_problem->getResidual() );

	// Neumann-Ulam solver for the correction.
 	d_nu_solver->linearProblem()->getLHS()->putScalar( 0.0 );
	d_nu_solver->walk();

	// Apply the correction.
	this->b_linear_problem->getLHS()->update( 
	    1.0, *(d_nu_solver->linearProblem()->getLHS()), 1.0 );

	// Update residual norm for convergence.
	this->b_linear_problem->computeResidual();
	this->b_linear_operator_split->applyInvM(
	    this->b_linear_problem->getResidual(),
	    this->b_linear_problem->getResidual() );
	residual_norm = this->b_linear_problem->getResidual()->normInf();

	// Update iteration count.
	++(this->b_num_iters);

	// Print iteration data.
	if ( this->b_linear_problem->getOperator()->getComm()->getRank() == 0 )
	{
	    std::cout << "Sequential MC Iteration " << this->b_num_iters 
		      << ": Residual = " << residual_norm << std::endl;
	}
	this->b_linear_problem->getOperator()->getComm()->barrier();
    }

    // Check for convergence.
    if ( residual_norm <= convergence_criteria )
    {
	this->b_is_converged = true;
    }
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//

#endif // end Chimera_SEQUENTIALMC_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_SequentialMC_def.hpp
//---------------------------------------------------------------------------//
