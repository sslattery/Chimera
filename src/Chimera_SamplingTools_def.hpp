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
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Ptr.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Stratify sample a global defined PDF to get the number of histories
 * required for each local state.
 */
template<class Scalar, class LO, class GO>
Teuchos::ArrayRCP<GO> SamplingTools::stratifySampleGlobalPDF( 
    const GO global_num_histories,
    const Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO> >& pdf )
{
    Teuchos::ArrayRCP<const Scalar> local_values = pdf->get1dView();
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType local_sum =
	std::accumulate( local_values.begin(), local_values.end(), 0.0 );

    typename Teuchos::ScalarTraits<Scalar>::magnitudeType global_sum =
	pdf->norm1();

    GO local_num_histories = std::floor( global_num_histories * 
					 local_sum / global_sum );

    // First, stratify sample the global pdf to get the number of histories
    // each local proc will generate.
    Teuchos::RCP<const Teuchos::Comm<int> > comm = pdf->getMap()->getComm();
    GO global_sum_check = 0;
    Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_num_histories, 
			Teuchos::Ptr<GO>(&global_sum_check) );
    
    GO global_remainder = global_num_histories - global_sum_check;
    if ( global_remainder > 0 )
    {
	if ( Teuchos::as<GO>(comm->getRank()) < global_remainder )
	{
	    ++local_num_histories;
	}
    }
    else if ( global_remainder < 0 )
    {
	if ( Teuchos::as<GO>(comm->getRank()) < std::abs(global_remainder) )
	{
	    --local_num_histories;
	}
    }
    comm->barrier();

    // Check that we maintained the global number of histories requested.
    remember( 
	global_sum_check = 0;
	Teuchos::reduceAll( *comm, Teuchos::REDUCE_SUM, local_num_histories, 
			    Teuchos::Ptr<GO>(&global_sum_check) );
	);

    testPostcondition( global_sum_check == global_num_histories );

    // Second, stratify sample the local pdf to get the number of histories to
    // be generated in each local state.
    typename Teuchos::ArrayRCP<GO> bin_histories( local_values.size() );
    typename Teuchos::ArrayRCP<Scalar>::const_iterator local_values_it;
    typename Teuchos::ArrayRCP<GO>::iterator bin_histories_it;
    for ( bin_histories_it = bin_histories.begin(),
	   local_values_it = local_values.begin();
	  bin_histories_it != bin_histories.end();
	  ++bin_histories_it, ++local_values_it )
    {
	*bin_histories_it = std::floor( local_num_histories * 
					(*local_values_it) / local_sum );
    }

    GO local_histories_sum = std::accumulate( bin_histories.begin(), 
					      bin_histories.end(), 0.0 );
    GO local_remainder = local_num_histories - local_histories_sum;
    GO local_index = 0;
    for ( bin_histories_it = bin_histories.begin();
	  bin_histories_it != bin_histories.end();
	  ++bin_histories_it )
    {
	local_index = std::distance( bin_histories.begin(), bin_histories_it );	

	if ( local_remainder > 0 )
	{
	    if (  local_index < local_remainder )
	    {
		++(*bin_histories_it);
	    }
	}
	else if ( local_remainder < 0 )
	{
	    if ( local_index < std::abs(local_remainder) )
	    {
		--(*bin_histories_it);
	    }
	}
    }

    // Check that we maintained the local number of histories requested.
    remember( 
	local_histories_sum = std::accumulate( bin_histories.begin(), 
					       bin_histories.end(), 0.0 );
	);

    testPostcondition( local_histories_sum == local_num_histories );
    
    // Return the number of histories to be generated in each local bin of the
    // PDF.
    return bin_histories;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Random sample a local discrete PDF for a new local state index.
 */
template<class Scalar, class LO, class RNG>
LO SamplingTools::sampleLocalDiscretePDF( 
    const Teuchos::ArrayView<Scalar>& pdf_values,
    const Teuchos::ArrayView<LO>& pdf_indices,
    const Teuchos::RCP<RNG>& rng )
{
    Scalar pdf_sum = std::accumulate( 
	pdf_values.begin(), pdf_values.end(), 0.0 );

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