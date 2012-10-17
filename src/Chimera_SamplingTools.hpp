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
 * \file Chimera_SamplingTools.hpp
 * \author Stuart R. Slattery
 * \brief Sampling tools declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_SAMPLINGTOOLS_HPP
#define Chimera_SAMPLINGTOOLS_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include <Tpetra_Vector.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class SamplingTools.
 * \brief A stateless class of tools for distribution function sampling.
 */
//---------------------------------------------------------------------------//
class SamplingTools
{
  public:

    //! Constructor.
    SamplingTools()
    { /* ... */ }

    //! Destructor.
    ~SamplingTools()
    { /* ... */ }

    // Stratify sample a global PDF.
    template<class Scalar, class LO, class GO>
    static Teuchos::ArrayRCP<GO> stratifySampleGlobalPDF( 
	const GO global_num_histories,
	const Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> >& pdf );

    // Random sample a local discrete PDF for a new state index.
    template<class Scalar, class LO, class RNG>
    static LO sampleLocalDiscretePDF( 
	const Teuchos::ArrayView<Scalar>& pdf_values,
	const Teuchos::ArrayView<LO>& pdf_indices,
	const Teuchos::RCP<RNG>& rng );
};

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_SamplingTools_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_SAMPLINGTOOLS_HPP

//---------------------------------------------------------------------------//
// end Chimera_SamplingTools.hpp
//---------------------------------------------------------------------------//
