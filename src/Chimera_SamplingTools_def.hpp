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
 * \file Chimera_SamplingTools_def.hpp
 * \author Stuart R. Slattery
 * \brief Sampling tools definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_SAMPLINGTOOLS_DEF_HPP
#define Chimera_SAMPLINGTOOLS_DEF_HPP

#include <numeric>

#include "Chimera_Assertion.hpp"
#include "Chimera_RNGTraits.hpp"

#include <Teuchos_as.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Stratify sample a global defined PDF to get the number of histories
 * required for each local state.
 */
template<class Scalar, class LO, class GO>
void SamplingTools::stratifySampleGlobalPDF( 
    const GO global_num_histories,
    const Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> >& pdf,
    Teuchos::ArrayRCP<LO>& local_histories_per_bin )
{
    Teuchos::ArrayRCP<Scalar> local_values = pdf->get1dView();
    typename Teuchos::ScalarTraits<Scalar>::magnitude_type local_sum =
	std::accumulate( local_values.begin(), local_values.end(), 0.0 );

    typename Teuchos::ScalarTraits<Scalar>::magnitude_type global_sum =
	pdf->norm1();

    GO local_num_histories = global_num_histories * 
			     std::floor( local_sum / global_sum );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Random sample a local discrete PDF for a new local state index.
 */
template<class Scalar, class LO, class RNG>
LO SamplingTools::sampleLocalDiscretePDF( 
    const Scalar pdf_sum,
    const Teuchos::ArrayView<Scalar>& pdf_values,
    const Teuchos::ArrayView<LO>& pdf_indices,
    const Teuchos::RCP<RNG>& rng )
{
    Scalar zeta = pdf_sum * 
		  ( Teuchos::as<Scalar>(RNGTraits<RNG>::generate(*rng)) /
		    Teuchos::as<Scalar>(RNGTraits<RNG>::max(*rng)) );

    Scalar cdf = 0.0;
    LO new_state_index = 0;
    typename Teuchos::ArrayView<Scalar>::const_iterator value_begin =
	pdf_values.begin();
    typename Teuchos::ArrayView<Scalar>::const_iterator value_iterator;
    for ( value_iterator = pdf_values.begin();
	  value_iterator != pdf_values.end();
	  ++value_iterator )
    {
	cdf += *value_iterator;
	if ( zeta <= cdf )
	{
	    new_state_index = Teuchos::as<LO>(
		std::distance( value_begin, value_iterator ) );

	    return pdf_indices[ new_state_index ];
	}
    }

    testPostcondition( zeta <= cdf );
    return 0;
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_SAMPLINGTOOLS_HPP

//---------------------------------------------------------------------------//
// end Chimera_SamplingTools.hpp
//---------------------------------------------------------------------------//
