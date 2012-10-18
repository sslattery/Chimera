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
 * \file Chimera_HistoryBuffer_def.hpp
 * \author Stuart R. Slattery
 * \brief HistoryBuffer definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_HISTORYBUFFER_DEF_HPP
#define Chimera_HISTORYBUFFER_DEF_HPP

#include "Chimera_Assertion.hpp"

#include <Teuchos_ArrayView.hpp>

#include <Tpetra_Distributor.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<class HT>
HistoryBuffer<HT>::HistoryBuffer( const RCP_TpetraMap& state_map )
    : d_state_map( state_map )
    , d_buffer( 0 )
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
template<class HT>
HistoryBuffer<HT>::~HistoryBuffer()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Communicate the buffer to its destinations and return a fresh
 * history bank.
*/
template<class HT>
HistoryBank<HT> HistoryBuffer<HT>::communicate()
{
    // Get the global states for the histories in the buffer.
    Teuchos::Array<GO> global_states( d_buffer.size() );
    typename Teuchos::Array<GO>::iterator global_state_it;
    typename Teuchos::Array<HT>::const_iterator buffer_it;
    for ( buffer_it = d_buffer.begin(), global_state_it = global_states.begin();
	  buffer_it != d_buffer.end(); 
	  ++buffer_it, ++global_state_it )
    {
	*global_state_it = buffer_it->globalState();
    }

    // Get the destination procs for those global states.
    Teuchos::Array<int> destination_procs( d_buffer.size() );
    Tpetra::LookupStatus lookup_status = 
	d_state_map->getRemoteIndexList( global_states(), destination_procs() );
    testInvariant( lookup_status == Tpetra::AllIDsPresent );

    global_states.clear();

    // Redistribute the histories to their destinations.
    Tpetra::Distributor distributor( d_state_map->getComm() );
    GO num_incoming_histories = 
	distributor.createFromSends( destination_procs() );

    destination_procs.clear();

    Teuchos::Array<HT> incoming_buffer( num_incoming_histories );
    Teuchos::ArrayView<const HT> outgoing_buffer = d_buffer();
    distributor.doPostsAndWaits( outgoing_buffer, 1, incoming_buffer() );

    // Return a new history bank populated with the incoming buffer.
    return HistoryBank<HT>( incoming_buffer );
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_HISTORYBUFFER_DEF_HPP

//---------------------------------------------------------------------------//
// end Chimera_HistoryBuffer_def.hpp
//---------------------------------------------------------------------------//
