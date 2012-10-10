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
 * \file Chimera_LinearOperatorSplit.hpp
 * \author Stuart Slattery
 * \brief Linear operator split interface definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_LINEAROPERATORSPLIT_HPP
#define Chimera_LINEAROPERATORSPLIT_HPP

#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*! 
 * \class LinearOperatorSplit
 * \brief Interface definition for linear operator splittings (A = M - N).
 */
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO>
class LinearOperatorSplit
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                  scalar_type;
    typedef LO                                      local_ordinal_type;
    typedef GO                                      global_ordinal_type;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO>         TpetraCrsMatrix;
    typedef Teuchos::RCP<TpetraCrsMatrix>           RCP_TpetraCrsMatrix;
    typedef Tpetra::Vector<Scalar,LO,GO>            TpetraVector;
    typedef Teuchos::RCP<TpetraVector>              RCP_TpetraVector;
    //@}

    //! Constructor.
    LinearOperatorSplit()
    { /* ... */ }

    //! Destructor.
    virtual ~LinearOperatorSplit()
    { /* ... */ }

    //! Split the operator.
    virtual void split() = 0;

    //! Get the linear operator (A).
    RCP_TpetraCrsMatrix linearOperator()
    { return b_linear_operator; }

    //! Get the iteration matrix (M^-1 N).
    RCP_TpetraCrsMatrix iterationMatrix()
    { return b_iteration_matrix; }

    //! Apply M^-1 to a vector (M^-1 x = y).
    virtual void 
    applyInvM( const RCP_TpetraVector& x, RCP_TpetraVector& y ) = 0;

  protected:

    // Original linear operator.
    RCP_TpetraCrsMatrix b_linear_operator;

    // Iteration matrix.
    RCP_TpetraCrsMatrix b_iteration_matrix;
};

} // end namespace Chimera

#endif // end Chimera_LINEAROPERATORSPLIT_HPP

//---------------------------------------------------------------------------//
// end Chimera_LinearOperatorSplit.hpp
//---------------------------------------------------------------------------//

