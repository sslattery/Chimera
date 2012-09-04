//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_BCStrategy_Factory.hpp
 * \author Stuart R. Slattery
 * \brief Boundary condition factory for transient poisson problem.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_BCSTRATEGY_FACTORY_HPP
#define TRANSIENTPOISSON_BCSTRATEGY_FACTORY_HPP

#include "TransientPoisson_BCStrategy_Dirichlet_Constant.hpp"

#include <Teuchos_RCP.hpp>

#include <Panzer_Traits.hpp>
#include <Panzer_BCStrategy_TemplateManager.hpp>
#include <Panzer_BCStrategy_Factory.hpp>
#include <Panzer_BCStrategy_Factory_Defines.hpp>

namespace Chimera
{
namespace TransientPoisson
{

PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER( "Constant",
					    BCStrategy_Dirichlet_Constant,
					    BCStrategy_Dirichlet_Constant )

struct BCStrategyFactory : public panzer::BCStrategyFactory
{

    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy( const panzer::BC& bc, 
		     const Teuchos::RCP<panzer::GlobalData>& global_data ) const
    {
	Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs_tm =
	    Teuchos::rcp( new panzer::BCStrategy_TemplateManager<panzer::Traits> );

	bool found = false;

	PANZER_BUILD_BCSTRATEGY_OBJECTS( 
	    "Constant",
	    Chimera::TransientPoisson::BCStrategy_Dirichlet_Constant,
	    BCStrategy_Dirichlet_Constant );

	TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, 
				   "Error - the BC Strategy called \"" << bc.strategy() <<
				   "\" is not a valid identifier in the BCStrategyFactory.  Either add a "
				   "valid implementation to your factory or fix your input file.  The "
				   "relevant boundary condition is:\n\n" << bc << std::endl);
      
	return bcs_tm;
    }
};

} // end namespace TransientPoisson
} // end namespace Chimera

#endif // end TRANSIENTPOISSON_BCSTRATEGY_FACTORY_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_BCStrategy_Factory.hpp
//---------------------------------------------------------------------------//


