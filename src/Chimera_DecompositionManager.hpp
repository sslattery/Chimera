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
// \file Chimera_DecompositionManager.hpp
// \author Stuart R. Slattery
// \brief Multiple-Set Overlapping-Domain Decomposition Manager declaration.
//---------------------------------------------------------------------------//

#ifndef Chimera_DECOMPOSITIONMANAGER_HPP
#define Chimera_DECOMPOSITIONMANAGER_HPP

#include "Chimera_LinearProblem.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Exporter.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class DecompositionManager
 * \brief Manager for Multiple-Set Overlapping-Domain Decomposition.
 */
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO>
class DecompositionManager
{
  public:

    //@{
    //! Typdefs
    typedef Scalar                                        scalar_type;
    typedef LO                                            local_ordinal_type;
    typedef GO                                            global_ordinal_type;
    typedef Teuchos::RCP<LinearProblem<Scalar,LO,GO> >    RCP_LinearProblem;
    typedef Teuchos::RCP<Teuchos::ParameterList>          RCP_ParameterList;
    typedef Teuchos::RCP<Tpetra::Export<LO,GO> >          RCP_TpetraExport;
    //@}

    // Constructor.
    DecompositionManager( const RCP_LinearProblem& base_problem,
			  const RCP_ParameterList& plist );

    // Deconstructor.
    ~DecompositionManager();

    // Export the base linear problem data to the MSOD decomposition.
    void exportBaseDataToDecomposition();

    // Export the decomposed linear problem data to the base decomposition.
    void exportDecomposedDataToBase();

    // Get the base linear problem.
    RCP_LinearProblem getBaseProblem()
    { return d_base_problem; }

    // Get the MSOD decomposed linear problem.
    RCP_LinearProblem getDecomposedProblem()
    { return d_msod_problem; }

    // Get the overlap fraction.
    double getOverlapFraction() const
    { return d_overlap_fraction; }

    // Get the number of sets.
    int getNumSets() const
    { return d_num_sets; }

  private:

    // Base linear problem.
    RCP_LinearProblem d_base_problem;

    // Domain overlap fraction.
    double d_overlap_fraction;

    // Number of sets.
    int d_num_sets;

    // MSOD decomposed problem.
    RCP_LinearProblem d_msod_problem;

    // Base-to-MSOD exporter.
    RCP_TpetraExport d_base_to_msod_export;

    // MSOD-to-Base exporter.
    RCP_TpetraExport d_msode_to_base_export;
};

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_DecompositionManager_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_DECOMPOSITIONMANAGER_HPP

//---------------------------------------------------------------------------//
// end Chimera_DecompositionManager.hpp
//---------------------------------------------------------------------------//

