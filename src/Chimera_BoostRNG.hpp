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
// \file Chimera_BoostRNG.hpp
// \author Stuart R. Slattery
// \brief RNGTraits implementation for Boost Psuedo Random Number Generators.
//---------------------------------------------------------------------------//

#ifndef Chimera_BOOSTRNG_HPP
#define Chimera_BOOSTRNG_HPP

#include "Chimera_RNGTraits.hpp"

namespace Chimera
{
namespace Solvers
{
//---------------------------------------------------------------------------//
template<>
template<class BoostRNG>
class RNGTraits<BoostRNG>
{
    typedef RNGType        BoostRNG;
    typedef result_type    BoostRNG::result_type;

    static inline result_type generate( const BoostRNG& boost_rng )
    { return boost_rng(); }

    static inline void reset( const BoostRNG& boost_rng )
    { boost_rng.reset(); }

    static inline void setSeed( const BoostRNG& boost_rng, const unsigned int seed )
    { boost_rng.seed( seed ); }

    static inline result_type min( const BoostRNG& boost_rng )
    { return boost_rng.min(); }

    static inline result_type max( const BoostRNG& boost_rng )
    { return boost_rng.max(); }
};

//---------------------------------------------------------------------------//

} // end namespace Solvers

} // end namespace Chimera

#endif // end Chimera_BOOSTRNG_HPP

//---------------------------------------------------------------------------//
// end Chimera_BoostRNG.hpp
//---------------------------------------------------------------------------//

