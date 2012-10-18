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
 * \file Chimera_HistoryBank.hpp
 * \author Stuart R. Slattery
 * \brief HistoryBank definition.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_BANK_HPP
#define Chimera_BANK_HPP

#include "Chimera_Assertion.hpp"

#include <Teuchos_Array.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class HistoryBank
 * \brief History bank. 
 *
 * This class is a stack using a Teuchos::Array as the underlying data
 * structure so that we may use it for communication operations.
 */
//---------------------------------------------------------------------------//
template<class HT>
class HistoryBank
{
  public:

    //@{
    //! Typedefs.
    typedef HT                                             history_type;
    typedef typename Teuchos::Array<HT>::size_type         size_type;
    //@}

    //! Constructor.
    HistoryBank()
    { /* ... */ }

    //! Destructor.
    ~HistoryBank()
    { /* ... */ }

    //! Set the history stack.
    void setStack( const Teuchos::Array<HT>& histories ) 
    { d_histories = history_stack; }

    //! Return if the bank is empty.
    bool empty () const
    { return d_histories.empty(); }

    //! Return the number of histories left in the bank.
    size_type size() const
    { return d_histories.size(); }

    //! Access the top history in the stack.
    inline const HT& top() const;

    //! Push history onto the stack.
    void push( const HT& history )
    { d_histories.push_back(history); }

    //! Pop a history off of the stack.
    inline HT pop();

  private:

    // History stack.
    Teuchos::Array<HT> d_histories;
};

//---------------------------------------------------------------------------//
// Inline functions.
//---------------------------------------------------------------------------//
template<class HT>
const HT& HistoryBank<HT>::top() const
{
    testPrecondition( !empty() );
    return d_histories.back();
}

//---------------------------------------------------------------------------//
template<class HT>
HT HistoryBank<HT>::pop()
{
    testPrecondition( !empty() );
    HT history = d_histories.back();
    d_histories.pop_back();
    return history;
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_BANK_HPP

//---------------------------------------------------------------------------//
// end Chimera_HistoryBank.hpp
//---------------------------------------------------------------------------//

