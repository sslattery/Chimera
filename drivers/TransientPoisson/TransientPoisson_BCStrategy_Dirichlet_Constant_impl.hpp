//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_BCStrategy_Dirichlet_Constant_impl.hpp
 * \author Stuart R. Slattery
 * \brief Constant Dirichlet boundary condition implementation.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_BCSTRATEGY_DIRICHLET_CONSTANT_IMPL_HPP
#define TRANSIENTPOISSON_BCSTRATEGY_DIRICHLET_CONSTANT_IMPL_HPP

#include <vector>

#include <Teuchos_Assert.hpp>

namespace Chimera
{
namespace TransientPoisson
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename EvaluationType>
BCStrategy_Dirichlet_Constant<EvaluationType>::BCStrategy_Dirichlet_Constant( 
    const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data )
    : panzer::BCStrategy_Dirichlet_DefaultImpl<EvaluationType>(bc,global_data)
{
    TEUCHOS_ASSERT( this->m_bc.strategy() == "Constant" );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Setup the base class data.
 */
template<typename EvaluationType>
void BCStrategy_Dirichlet_Constant<EvaluationType>::setup( 
    const panzer::PhysicsBlock& side_pb, 
    const Teuchos::ParameterList& user_data )
{
    // Get the dofs for the residual.
    this->required_dof_names.push_back( this->m_bc.equationSetName() );

    // Get the name of the residual.
    this->residual_name = "Residual_" + this->m_bc.identifier();

    // Map the residual to the dofs.
    this->residual_to_dof_names_map[residual_name] = this->m_bc.equationSetName();

    // Map the residual to the target.
    this->residual_to_target_field_map[residual_name] = 
	"Constant_" + this->m_bc.equationSetName();

    // Get the basis for the dofs.
    const std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >&
	dofs = side_pb.getProvidedDOFs();
    std::vector<std::pair<std::string,Teuchos::RCP<panzer::PureBasis> > >::const_iterator
	dof_iterator;
    for ( dof_iterator = dofs.begin(); dof_iterator != dofs.end(); ++dof_iterator )
    {
	if ( dof_iterator->first == this->m_bc.equationSetName() )
	{
	    this->basis = dof_iterator->second;
	}
    }

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis), std::runtime_error,
			       "Error the name \"" << this->m_bc.equationSetName()
			       << "\" is not a valid DOF for the boundary condition:\n"
			       << this->m_bc << "\n");
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build boundary condition evaluators.
 */
template<typename EvaluationType>
void BCStrategy_Dirichlet_Constant<EvaluationType>::buildAndRegisterEvaluators( 
    PHX::FieldManager<panzer::Traits>& fm,
    const panzer::PhysicsBlock& pb,
    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
    const Teuchos::ParameterList& models,
    const Teuchos::ParameterList& user_data ) const
{
    // Build a constant value evaluator.
    {
	Teuchos::ParameterList plist( "BC Constant Dirichlet" );
	plist.set( "Name", "Constant_" + this->m_bc.equationSetName() );
	plist.set( "Data Layout", basis->functional );
	plist.set( "Value", this->m_bc.params()->get<double>( "Value" ) );

	Teuchos::RCP<PHX::Evaluator<panzer::Traits> > dirichlet_op =
	    Teuchos::rcp(new panzer::Constant<EvaluationType,panzer::Traits>(plist) );

	fm.registerEvaluator<EvaluationType>( dirichlet_op );
    }
}

//---------------------------------------------------------------------------//

} // end namespace TransientPoisson
} // end namespace Chimera

//---------------------------------------------------------------------------//

#endif // end TRANSIENTPOISSON_BCSTRATEGY_DIRICHLET_CONSTANT_IMPL_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_BCStrategy_Dirichlet_Constant_impl.hpp
//---------------------------------------------------------------------------//

