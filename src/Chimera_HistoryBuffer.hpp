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
 * \file Chimera_HistoryBuffer.hpp
 * \author Stuart R. Slattery
 * \brief HistoryBuffer declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_HISTORYBUFFER_HPP
#define Chimera_HISTORYBUFFER_HPP

#include "Chimera_HistoryBank.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_Map.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class HistoryBuffer
 * \brief History buffer for communicating histories.
 */
//---------------------------------------------------------------------------//
template<class HT>
class HistoryBuffer
{
  public:

    //@{
    //! Typedefs.
    typedef HT                                                history_type;
    typedef typename Teuchos::Array<HT>::size_type            size_type;
    //@}

    //! Constructor.
    HistoryBuffer();

    //! Destructor.
    ~HistoryBuffer();

    //! Push back a history into the buffer.
    void pushBack( const HT& outgoing_history )
    { d_buffer.push_back( outgoing_history ); }

    // Communicate the history buffer to its destinations and return a fresh
    // bank of local histories.
    template<class LO, class GO>
    HistoryBank<HT> communicate( 
	const Teuchos::RCP<const Tpetra::Map<LO,GO> >& state_map );

    //! Flush the buffer.
    void flush()
    { d_buffer.clear(); }

    //! Return if the buffer is empty.
    bool empty() const
    { return d_buffer.empty(); }

    //! Return the number of histories in the buffer.
    size_type size() const
    { return d_buffer.size(); }

  private:

    // History buffer.
    Teuchos::Array<HT> d_buffer;
};

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_HistoryBuffer_def.hpp"

//---------------------------------------------------------------------------//

#endif // end Chimera_HISTORYBUFFER_HPP

//---------------------------------------------------------------------------//
// end Chimera_HistoryBuffer.hpp
//---------------------------------------------------------------------------//

