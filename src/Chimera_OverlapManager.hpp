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
// \file Chimera_OverlapManager.hpp
// \author Stuart R. Slattery
// \brief Multiple-Set Overlapping-Domain Decomposition Manager declaration.
//---------------------------------------------------------------------------//

#ifndef Chimera_OVERLAPMANAGER_HPP
#define Chimera_OVERLAPMANAGER_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Export.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
  \class OverlapManager
  \brief Manager for Overlapping-Domain Decomposition in Neumann-Ulam solver.
*/
//---------------------------------------------------------------------------//
template<class Scalar, class LO, class GO>
class OverlapManager
{
  public:

    //@{
    //! Typdefs
    typedef Scalar                                        scalar_type;
    typedef LO                                            local_ordinal_type;
    typedef GO                                            global_ordinal_type;
    typedef Teuchos::RCP<Teuchos::ParameterList>          RCP_ParameterList;
    typedef Tpetra::Vector<Scalar,LO,GO>                  TpetraVector;
    typedef Teuchos::RCP<TpetraVector>                    RCP_TpetraVector;
    typedef Tpetra::CrsMatrix<Scalar,LO,GO>               TpetraCrsMatrix;
    typedef Teuchos::RCP<TpetraCrsMatrix>                 RCP_TpetraCrsMatrix;
    typedef Teuchos::RCP<Tpetra::Export<LO,GO> >          RCP_TpetraExport;
    //@}

    // Constructor.
    OverlapManager( const RCP_TpetraCrsMatrix& iteration_matrix,
		    const RCP_TpetraCrsMatrix& probability_matrix,
		    const RCP_TpetraVector& lhs,
		    const RCP_ParameterList& plist );

    // Deconstructor.
    ~OverlapManager();

    // Determine if a global state is in the overlap owned by this proc.
    bool isOverlapGlobalElement( const GO global_state );

    // Get the number of overlapping states.
    GO getNumOverlap() const
    { return d_num_overlap; }

    // Get the overlapping iteration matrix.
    RCP_TpetraCrsMatrix getOverlapIterationMatrix() 
    { return d_overlap_iteration_matrix; }

    // Get the overlapping probability matrix.
    RCP_TpetraCrsMatrix getOverlapProbabilityMatrix() 
    { return d_overlap_probability_matrix; }

    // Get the overlapping LHS.
    RCP_TpetraVector getOverlappingLHS()
    { return d_overlapping_lhs; }

    // Get the overlap-to-base export.
    RCP_TpetraExport getOverlapToBaseExport()
    { return d_overlap_to_base_export; }

  private:

    // Build the overlap.
    void buildOverlap();

  private:

    // Overlapping iteration matrix.
    RCP_TpetraCrsMatrix d_overlap_iteration_matrix;

    // Overlapping probability matrix.
    RCP_TpetraCrsMatrix d_overlap_probability_matrix;

    // Overlapping LHS.
    RCP_TpetraVector d_overlap_lhs;

    // Number of overlapping states.
    GO d_num_overlap;

    // Overlap-to-Base exporter.
    RCP_TpetraExport d_overlap_to_base_export;
};

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_OverlapManager_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_OVERLAPMANAGER_HPP

//---------------------------------------------------------------------------//
// end Chimera_OverlapManager.hpp
//---------------------------------------------------------------------------//

