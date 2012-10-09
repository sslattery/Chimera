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
 * \file Chimera_JacobiSplit.hpp
 * \author Stuart R. Slattery
 * \brief Jacobi split declaration.
 */ 
//---------------------------------------------------------------------------//

#ifndef Chimera_JACOBISPLIT_HPP
#define Chimera_JACOBISPLIT_HPP

#include "Chimera_LinearOperatorSplit.hpp"

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class JacobiSplit
 * \brief Jacobi matrix splitting implementation. 
 *
 * This class implements either the Jacobi stationary method or a Jacobi
 * preconditioned Richardson iteration.
 */
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO>
class JacobiSplit : public LinearOperatorSplit<Scalar,LO,GO>
{
  public:

    //@{
    //! Typdefs
    typedef Scalar                                  scalar_type;
    typedef LO                                      local_ordinal_type;
    typedef GO                                      global_ordinal_type;
    typedef LinearOperatorSplit<Scalar,LO,GO>       Base;
    typedef Base::RCP_TpetraCrsMatrix               RCP_TpetraCrsMatrix;
    typedef Base::RCP_TpetraVector                  RCP_TpetraVector;
    //@}

    // Constructor.
    JacobiSplit( const RCP_TpetraCrsMatrix& linear_operator );

    // Destructor.
    ~JacobiSplit();

    // Split the operator with Jacobi splitting.
    void split();

    // Apply M^-1 to a vector (M^-1 x = y).
    void applyInvM( const RCP_TpetraVector& x, RCP_TpetraVector& y );
};

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_JacobiSplit_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_JACOBISPLIT_HPP

//---------------------------------------------------------------------------//
// end Chimera_JacobiSplit.hpp
//---------------------------------------------------------------------------//

