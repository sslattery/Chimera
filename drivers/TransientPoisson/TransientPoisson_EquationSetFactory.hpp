//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_EquationSetFactory.hpp
 * \author Stuart R. Slattery
 * \brief Equation set factory for the transient Poisson problem.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_EQSET_FACTORY_HPP
#define TRANSIENTPOISSON_EQSET_FACTORY_HPP

#include "TransientPoisson_EquationSet.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_TestForException.hpp>

#include <Panzer_EquationSet_Factory.hpp>
#include <Panzer_EquationSet_Factory_Defines.hpp>
#include <Panzer_InputEquationSet.hpp>
#include <Panzer_CellData.hpp>
#include <Panzer_GlobalData.hpp>

//---------------------------------------------------------------------------//

namespace Chimera
{
namespace TransientPoissson
{

// Define a class for building the equation sets with a macro.
PANZER_DECLARE_EQSET_TEMPLATE_BUILDER( "TransientPoisson", 
				       TransientPoissonEquationSet, 
				       TransientPoissonEquationSet );

// Define the factory method for the equations.
class EquationSetFactory : public panzer::EquationSetFactory
{
  public:

    // Factory Method. Note here that the variables must be named this
    // way. The panzer macros for building objects expect these.
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet( const panzer::InputEquationSet& ies,
		      const panzer::CellData& cell_data,
		      const Teuchos::RCP<panzer::GlobalData>& global_data,
		      const bool build_transient_support ) const
    {
	// Setup a template manager.
	Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
	    eq_set = Teuchos::rcp(
		new panzer::EquationSet_TemplateManager<panzer::Traits> );

	// Build the equation set.
	bool found = false;
	PANZER_BUILD_EQSET_OBJECTS( "TransientPoisson",
				    TransientPoissonEquationSet, 
				    TransientPoissonEquationSet );

	// Verify that it was valid.
	TEUCHOS_TEST_FOR_EXCEPTION( 
	    found, std::logic_error,
	    "Failed to find TransientPoisson equation set" );

	return eq_set;
    }
};

} // end namespace TransientPoisson
} // end namespace Chimera

#endif // end TRANSIENTPOISSON_EQSET_FACTORY_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_EquationSetFactory.hpp
//---------------------------------------------------------------------------//

