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
 * \file Chimera_LinearProblem.hpp
 * \author Stuart Slattery
 * \brief Linear problem declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_LINEARPROBLEM_HPP
#define Chimera_LINEARPROBLEM_HPP

#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
// \class LinearProblem
// \brief Linear problem container for Chimera.
//---------------------------------------------------------------------------//
template<class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal>
class LinearProblem
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                               scalar_type;
    typedef LocalOrdinal                                         local_ordinal_type;
    typedef GlobalOrdinal                                        global_ordinal_type;
    typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal>    TpetraVector;
    typedef Teuchos::RCP<TpetraVector>                           RCP_TpetraVector;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> TpetraCrsMatrix;
    typedef Teuchos::RCP<TpetraCrsMatrix>                        RCP_TpetraCrsMatrix;
    //@}

  public:

    // Constructor.
    LinearProblem( const RCP_TpetraCrsMatrix& A, const RCP_TpetraVector& x, 
		   const RCP_TpetraVector& b );

    // Destructor.
    ~LinearProblem();

    //! Set the operator of the linear problem.
    void setOperator( const RCP_TpetraCrsMatrix& A )
    { d_A = A; }

    //! Set the solution vector of the linear problem.
    void setLHS( const RCP_TpetraVector& x )
    { d_x = x; }

    //! Set the right-hand sideof the linear problem.
    void setRHS( const RCP_TpetraVector& b )
    { d_b = b; }

    //! Get the operator of the linear problem.
    RCP_TpetraCrsMatrix getOperator()
    { return d_A; }

    //! Get the solution vector of the linear problem.
    RCP_TpetraVector getLHS()
    { return d_x; }

    //! Get the right-hand side of the linear problem.
    RCP_TpetraVector getRHS()
    { return d_b; }

    // Compute the residual of the linear problem.
    void computeResidual();

    //! Get the residual of the linear problem.
    RCP_TpetraVector getResidual()
    { return d_r; }

  private:
    
    // Linear operator.
    RCP_TpetraCrsMatrix d_A;

    // Solution vector.
    RCP_TpetraVector d_x;

    // Right-hand side.
    RCP_TpetraVector d_b;

    // Linear system residual.
    RCP_TpetraVector d_r;
};

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_LinearProblem_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_LINEARPROBLEM_HPP

//---------------------------------------------------------------------------//
// end Chimera_LinearProblem.hpp
//---------------------------------------------------------------------------//

