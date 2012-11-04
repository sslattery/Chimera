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
 * \file Chimera_StationarySolver_def.hpp
 * \author Stuart R. Slattery
 * \brief StationarySolver definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_STATIONARYSOLVER_DEF_HPP
#define Chimera_STATIONARYSOLVER_DEF_HPP

#include "Chimera_Assertion.hpp"
#include "Chimera_LinearOperatorSplitFactory.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_as.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO>
StationarySolver<Scalar,LO,GO>::StationarySolver( 
    const RCP_LinearProblem& linear_problem,
    const RCP_ParameterList& plist )
{
    testPrecondition( !linear_problem.is_null() );
    testPrecondition( !plist.is_null() );

    this->b_linear_problem = linear_problem;
    this->b_linear_operator_split = LinearOperatorSplitFactory::create( 
	plist, linear_problem->getOperator() );
    this->b_tolerance = plist->get<Scalar>("TOLERANCE");
    this->b_max_num_iters = plist->get<int>("MAX ITERS");
    this->b_num_iters = 0;
    this->b_is_converged = false;

    d_stationary_iteration = 
	Teuchos::rcp( new StationaryIteration<Scalar,LO,GO>( 
			  this->b_linear_problem,
			  this->b_linear_operator_split ) );

    testPostcondition( !this->b_linear_operator_split.is_null() );
    testPostcondition( !d_stationary_iteration.is_null() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Scalar, class LO, class GO>
StationarySolver<Scalar,LO,GO>::~StationarySolver()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Iterate until convergence.
 */
template<class Scalar, class LO, class GO>
void StationarySolver<Scalar,LO,GO>::iterate()
{
    // Generate the convergence criteria.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType b_norm = 
	this->b_linear_problem->getRHS()->normInf();
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
	convergence_criteria = this->b_tolerance * b_norm;

    // Setup.
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
	residual_norm = 1.0;
    this->b_num_iters = 0;
    this->b_is_converged = false;

    // Iterate.
    while ( residual_norm > 
	    Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(
		convergence_criteria) &&
	    this->b_num_iters < this->b_max_num_iters )
    {
	// Do a stationary iteration.
	d_stationary_iteration->doOneIteration();

	// Compute the new residual.
	this->b_linear_problem->computeResidual();

	this->b_linear_operator_split->applyInvM( 
	    this->b_linear_problem->getResidual(), 
	    this->b_linear_problem->getResidual() );

	residual_norm = this->b_linear_problem->getResidual()->normInf();

	++(this->b_num_iters);
    }

    // Check for convergence.
    if ( residual_norm <= 
	 Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(
	     convergence_criteria) )
    {
	this->b_is_converged = true;
    }
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_StationarySolver_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_STATIONARYSOLVER_HPP

//---------------------------------------------------------------------------//
// end Chimera_StationarySolver.hpp
//---------------------------------------------------------------------------//
