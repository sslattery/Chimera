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
// \file Chimera_RNGTraits.hpp
// \author Stuart R. Slattery
// \brief Traits definition for random number generators.
//---------------------------------------------------------------------------//

#ifndef Chimera_RNGTRAITS_HPP
#define Chimera_RNGTRAITS_HPP

#include <Teuchos_RCP.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Dummy struct. If a type does not create a specialization this will
 * not compile.
 */
template<typename UndefinedRNGType>
struct UndefinedRNGTraits
{
    static inline UndefinedRNGType notDefined() 
    { return UndefinedRNGType::this_type_is_missing_a_specialization(); }
};

//---------------------------------------------------------------------------//
/*!
  \class RNGTraits
  \brief Random number generator traits definitions.

  RNGTraits provide access to random number generators.
*/
//---------------------------------------------------------------------------//
template<class RNGType>
class RNGTraits
{
  public:

    //@{
    //! Typedef for rng type.
    typedef RNGType rng_type;

    //! Typedef for result type.
    typedef typename RNGType::result_type result_type;
    //@}

    /*!
     * \brief Create a random number generator in the base state.
     */
    static inline Teuchos::RCP<RNGType> create()
    { UndefinedRNGTraits<RNGType>::notDefined(); return 0; }

    /*! 
     * \brief Generate a random number.
     */
    static inline result_type generate( RNGType& rng )
    { UndefinedRNGTraits<RNGType>::notDefined(); return 0; }

    /*!
     * \brief Set the current seed of the random number generator.
     */
    static inline void setSeed( RNGType& rng, const result_type seed )
    { UndefinedRNGTraits<RNGType>::notDefined(); }

    /*!
     * \brief Lower bound of random range.
     */
    static inline result_type min( RNGType& rng )
    { UndefinedRNGTraits<RNGType>::notDefined(); return 0; }

    /*!
     * \brief Upper bound of random range.
     */
    static inline result_type max( RNGType& rng )
    { UndefinedRNGTraits<RNGType>::notDefined(); return 0; }
};

//---------------------------------------------------------------------------//

} // end namespace Chimera

#endif // end Chimera_RNGTRAITS_HPP

//---------------------------------------------------------------------------//
// end Chimera_RNGTraits.hpp
//---------------------------------------------------------------------------//

