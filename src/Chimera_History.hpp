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
template<class Scalar, class LO, class GO>
class History
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                    scalar_type;
    typedef LO                                        local_ordinal_type;
    typedef GO                                        global_ordinal_type;
    //@}

    // Defaut constructor.
    History();

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

    //! Set the local history state.
    void setLocalState( const LO local_state )
    { d_local_state = local_state; }

    //! Set the global history state.
    void setGlobalState( const GO global_state )
    { d_global_state = global_state; }

    //! Terminate the history.
    void terminate()
    { d_active = false; }

    //! Get the history weight.
    Scalar weight() const
    { return d_weight; }

    //! Get the local history state.
    LO localState() const
    { return d_local_state; }

    //! Get the global history state.
    GO globalState() const 
    { return d_global_state; }

    //! Get the active history state.
    bool active() const
    { return d_active; }

  private:

    // History weight.
    Scalar d_weight;

    // Local history state.
    LO d_local_state;

    // Global history state.
    GO d_global_state;

    // Active history state.
    bool d_active;
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
typedef Chimera::History<float,int,int>                         flt_int_int;
typedef Chimera::History<double,int,int>                        dbl_int_int;
typedef Chimera::History<float,int,long int>                    flt_int_long;
typedef Chimera::History<double,int,long int>                   dbl_int_long;
typedef Chimera::History<float,int,long long int>               flt_int_ll;
typedef Chimera::History<double,int,long long int>              dbl_int_ll;
typedef Chimera::History<float,long int,long int>               flt_long_long;
typedef Chimera::History<double,long int,long int>              dbl_long_long;
typedef Chimera::History<float,long int,long long int>          flt_long_ll;
typedef Chimera::History<double,long int,long long int>         dbl_long_ll;
typedef Chimera::History<float,long long int,long long int>     flt_ll_ll;
typedef Chimera::History<double,long long int,long long int>    dbl_ll_ll;

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_int_int>
    : public DirectSerializationTraits<Ordinal,flt_int_int >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_int_int>
    : public DirectSerializationTraits<Ordinal,dbl_int_int >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_int_long>
    : public DirectSerializationTraits<Ordinal,flt_int_long >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_int_long>
    : public DirectSerializationTraits<Ordinal,dbl_int_long >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_int_ll>
    : public DirectSerializationTraits<Ordinal,flt_int_ll >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_int_ll>
    : public DirectSerializationTraits<Ordinal,dbl_int_ll >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_long_long>
    : public DirectSerializationTraits<Ordinal,flt_long_long >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_long_long>
    : public DirectSerializationTraits<Ordinal,dbl_long_long >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_long_ll>
    : public DirectSerializationTraits<Ordinal,flt_long_ll >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_long_ll>
    : public DirectSerializationTraits<Ordinal,dbl_long_ll >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,flt_ll_ll>
    : public DirectSerializationTraits<Ordinal,flt_ll_ll >
{};

template<typename Ordinal>
class SerializationTraits<Ordinal,dbl_ll_ll>
    : public DirectSerializationTraits<Ordinal,dbl_ll_ll >
{};

} // end namespace Teuchos

//---------------------------------------------------------------------------//

#endif // end Chimera_HISTORY_HPP

//---------------------------------------------------------------------------//
// end Chimera_History.hpp
//---------------------------------------------------------------------------//

