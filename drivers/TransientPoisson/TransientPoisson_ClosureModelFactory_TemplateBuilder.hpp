//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_ClosureModelFactory_TemplateBuilder.hpp
 * \author Stuart R. Slattery
 * \brief ClosureModelFactory builder.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_CLOSUREMODEL_BUILDER_HPP
#define TRANSIENTPOISSON_CLOSUREMODEL_BUILDER_HPP

#include "TransientPoisson_ClosureModelFactory.hpp"

#include <Teuchos_RCP.hpp>

#include <Panzer_Base.hpp>

namespace Chimera
{
namespace TransientPoisson
{

class ClosureModelFactory_TemplateBuilder
{
  public:

    // Builder method.
    template<typename EvaluationType>
    Teuchos::RCP<panzer::Base> build() const
    {
	return Teuchos::rcp( static_cast<panzer::Base*>( 
				 new ClosureModelFactory<EvaluationType> ) );
    }

};

} // end namespace TransientPoisson
} // end namespace Chimera

#endif // end TRANSIENTPOISSON_CLOSUREMODEL_BUILDER_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_ClosureModelFactory_TemplateBuilder.hpp
//---------------------------------------------------------------------------//

