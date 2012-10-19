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
 * \file Chimera_History.hpp
 * \author Stuart R. Slattery
 * \brief History class declaration.
 */
//---------------------------------------------------------------------------//

#ifndef Chimera_HISTORY_HPP
#define Chimera_HISTORY_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_SerializationTraits.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \class History
 * \brief Encapsulation of a random walk history's state.
 */
//---------------------------------------------------------------------------//
template<class Scalar, class GO>
class History
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                    scalar_type;
    typedef GO                                        global_ordinal_type;
    //@}

    // Defaut constructor.
    History();

    // State constructor.
    History( Scalar weight, GO global_state );

    // Destructor.
    ~History();

    //! Set the history weight.
    void setWeight( const Scalar weight )
    { d_weight = weight; }

    //! Add to the history weight.
    void addWeight( const Scalar weight )
    { d_weight += weight; }

    //! Multiply the history weight.
    void multiplyWeight( const Scalar weight )
    { d_weight *= weight; }

    //! Set the global history state.
    void setGlobalState( const GO global_state )
    { d_global_state = global_state; }

    //! Get the history weight.
    Scalar weight() const
    { return d_weight; }

    //! Get the global history state.
    GO globalState() const 
    { return d_global_state; }

  private:

    // History weight.
    Scalar d_weight;

    // Global history state.
    GO d_global_state;
};

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "Chimera_History_def.hpp"

//---------------------------------------------------------------------------//
// SerializationTraits specializations.
//---------------------------------------------------------------------------//
namespace Teuchos
{
typedef Chimera::History<float,int>                         flt_int;
typedef Chimera::History<double,int>                        dbl_int;
typedef Chimera::History<float,long int>                    flt_long;
typedef Chimera::History<double,long int>                   dbl_long;
typedef Chimera::History<float,long long int>               flt_ll;
typedef Chimera::History<double,long long int>              dbl_ll;

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_int>
    : public DirectSerializationTraits<Ordinal,flt_int >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_int>
    : public DirectSerializationTraits<Ordinal,dbl_int >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_long>
    : public DirectSerializationTraits<Ordinal,flt_long >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_long>
    : public DirectSerializationTraits<Ordinal,dbl_long >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_ll>
    : public DirectSerializationTraits<Ordinal,flt_ll >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_ll>
    : public DirectSerializationTraits<Ordinal,dbl_ll >
{};

} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end Chimera_HISTORY_HPP

//---------------------------------------------------------------------------//
// end Chimera_History.hpp
//---------------------------------------------------------------------------//

