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
 * \file Chimera_StationaryIteration_def.hpp
 * \author Stuart R. Slattery
 * \brief StationaryIteration definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_STATIONARYITERATION_DEF_HPP
#define Chimera_STATIONARYITERATION_DEF_HPP

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class Scalar, class LO, class GO>
StationaryIteration<Scalar,LO,GO>::StationaryIteration(
    const RCP_LinearProblem& linear_problem,
    const RCP_LinearOperatorSplit& linear_operator_split )
    : d_linear_problem( linear_problem )
    , d_linear_operator_split( linear_operator_split )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Scalar, class LO, class GO>
StationaryIteration<Scalar,LO,GO>::~StationaryIteration()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Do one stationary iteration.
 */
template<class Scalar, class LO, class GO>
void StationaryIteration<Scalar,LO,GO>::doOneIteration()
{
    d_linear_operator_split->iterationMatrix()->apply( 
	*(d_linear_problem->getLHS()),
	*(d_linear_problem->getLHS()) );

    RCP_TpetraVector m_inv_b = Tpetra::createVector<Scalar,LO,GO>( 
	d_linear_problem->getOperator()->getRowMap() );

    d_linear_operator_split->applyInvM( d_linear_problem->getRHS(), m_inv_b );

    d_linear_problem->getLHS()->update( 1.0, *(m_inv_b), 1.0 );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_STATIONARYITERATION_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_StationaryIteration_def.hpp
//---------------------------------------------------------------------------//


