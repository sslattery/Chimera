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
 * \file Chimera_LinearProblem_def.hpp
 * \author Stuart Slattery
 * \brief Linear problem definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_LINEARPROBLEM_DEF_HPP
#define Chimera_LINEARPROBLEM_DEF_HPP

#include "Chimera_Assertion.hpp"

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for linear problem A*x = b.
 * \param A The linear operator.
 * \param x The solution vector.
 * \param b The right-hand side.
 */
template<class Scalar, class LO, class GO>
LinearProblem<Scalar,LO,GO>::LinearProblem( 
    const RCP_TpetraCrsMatrix& A, 
    const RCP_TpetraVector& x, 
    const RCP_TpetraVector& b )
    : d_A( A )
    , d_x( x )
    , d_b( b )
    , d_r( Tpetra::createVector<double,int>( d_A->getRowMap() ) )
{
    testPostcondition( !d_A.is_null() );
    testPostcondition( !d_x.is_null() );
    testPostcondition( !d_b.is_null() );
    testPostcondition( !d_r.is_null() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class Scalar, class LO, class GO>
LinearProblem<Scalar,LO,GO>::~LinearProblem()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Compute the residual of the linear problem r = b - A*x.
 */
template<class Scalar, class LO, class GO>
void LinearProblem<Scalar,LO,GO>::computeResidual()
{
    d_A->apply( *d_x, *d_r );
    d_r->update( 1.0, *d_b, -1.0 );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_LINEARPROBLEM_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_LinearProblem.hpp
//---------------------------------------------------------------------------//
