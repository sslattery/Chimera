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
 * \file Chimera_StationarySolver.hpp
 * \author Stuart R. Slattery
 * \brief StationarySolver declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_STATIONARYSOLVER_HPP
#define Chimera_STATIONARYSOLVER_HPP

#include <Chimera_LinearProblem.hpp>
#include <Chimera_StationaryIteration.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class StationarySolver
 * \brief Stationary method solver.
 */
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO>
class StationarySolver : public LinearSolver<Scalar,LO,GO>
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                    scalar_type;
    typedef LO                                        local_ordinal_type;
    typedef GO                                        global_ordinal_type;
    typedef LinearSolver<Scalar,LO,GO>                Base;
    typedef typename Base::RCP_LinearProblem          RCP_LinearProblem;
    typedef StationaryIteration<Scalar,LO,GO>         StationaryIterationType;
    typedef Teuchos::RCP<StationaryIterationType>     RCP_StationaryIteration;
    typedef Teuchos::RCP<Teuchos::ParameterList>      RCP_ParameterList;
    //@}

    //! Constructor.
    StationarySolver( const RCP_LinearProblem& linear_problem,
		      const RCP_ParameterList& plist );

    //! Destructor.
    ~StationarySolver();

    //! Iterate until convergence.
    void iterate();

  private:

    // Stationary iteration.
    RCP_StationaryIteration d_stationary_iteration;
};

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
