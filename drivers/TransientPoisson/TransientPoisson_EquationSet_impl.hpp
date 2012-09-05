//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_EquationSet_impl.hpp
 * \author Stuart R. Slattery
 * \brief Equation set implementation for transient Poisson problem.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_EQSET_IMPL_HPP
#define TRANSIENTPOISSON_EQSET_IMPL_HPP

#include "Teuchos_Assert.hpp"

#include <Panzer_Integrator_BasisTimesScalar.hpp>
#include <Panzer_Integrator_GradBasisDotVector.hpp>
#include <Panzer_Sum.hpp>

namespace Chimera
{
namespace TransientPoisson
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
template<typename EvaluationType>
TransientPoissonEquationSet<EvaluationType>::TransientPoissonEquationSet( 
    const panzer::InputEquationSet& ies,
    const panzer::CellData& cell_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    const bool build_transient_support )
:  panzer::EquationSet_DefaultImpl<EvaluationType>(
    ies, cell_data, global_data, build_transient_support )
{
    // Verify that there is no prefix and then set the prefix.
    TEUCHOS_ASSERT( ies.prefix == "" );
    this->m_eqset_prefix = "";

    // Declare temperature variables.
    this->m_dof_names->push_back("TEMPERATURE");
    this->m_dof_gradient_names->push_back("GRAD_TEMPERATURE");
    this->m_residual_names->push_back("RESIDUAL_TEMPERATURE");
    this->m_scatter_name = "Scatter_RESIDUAL_TEMPERATURE";

    // Set time derivative if transient support enabled.
    if ( this->m_build_transient_support )
    {
	this-> m_dof_time_derivative_names->push_back("DOT_TEMPERATURE");
    }

    // Build the basis functions and integration rules.
    this->setupDOFs( cell_data.baseCellDimension() );
}

//---------------------------------------------------------------------------//
/*!
 * \brief Build and register the evaluators with the field manager by
 * assembling the Poisson equation in weak form.
 */
template<typename EvaluationType>
void TransientPoissonEquationSet<EvaluationType>::
buildAndRegisterEquationSetEvaluators(
    PHX::FieldManager<panzer::Traits>& field_manager,
    const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
    const Teuchos::ParameterList& user_data ) const
{
    // Build the transient operator if transient support enabled.
    if ( this->m_build_transient_support )
    {
	Teuchos::ParameterList plist( "Transient Residual" );
	plist.set( "Residual Name", "RESIDUAL_TEMPERATURE_TRANSIENT_OP" );
	plist.set( "Value Name", "DOT_TEMPERATURE" );
	plist.set( "Basis", this->m_basis );
	plist.set( "IR", this->m_int_rule );
	plist.set( "Multiplier", 1.0 );

	Teuchos::RCP< PHX::Evaluator<panzer::Traits> > transient_op = 
	    Teuchos::rcp( 
		new panzer::Integrator_BasisTimesScalar<EvaluationType,panzer::Traits>( plist ) );

	field_manager.registerEvaluator<EvaluationType>( transient_op );
    }
 
    // Build the diffusion operator.
    {
	double thermal_conductivity = 
	    user_data.get<double>( "Thermal Conductivity" );

	Teuchos::ParameterList plist( "Diffusion Residual" );
	plist.set( "Residual Name", "RESIDUAL_TEMPERATURE_DIFFUSION_OP" );
	plist.set( "Flux Name", "GRAD_TEMPERATURE" );
	plist.set( "Basis", this->m_basis );
	plist.set( "IR", this->m_int_rule );
	plist.set( "Multiplier", thermal_conductivity );

	Teuchos::RCP< PHX::Evaluator<panzer::Traits> > diffusion_op =
	    Teuchos::rcp(
		new panzer::Integrator_GradBasisDotVector<EvaluationType,panzer::Traits>( plist ) );

	field_manager.registerEvaluator<EvaluationType>( diffusion_op );
    }

    // Build the source operator. The actual value of the source is generated
    // by the closure model.
    {
	Teuchos::ParameterList plist( "Source Residual" );
	plist.set( "Residual Name", "RESIDUAL_TEMPERATURE_SOURCE_OP" );
	plist.set( "Value Name", "SOURCE_TEMPERATURE" );
	plist.set( "Basis", this->m_basis );
	plist.set( "IR", this->m_int_rule );
	plist.set( "Multiplier", -1.0 );

	Teuchos::RCP< PHX::Evaluator<panzer::Traits> > source_op =
	    Teuchos::rcp(
		new panzer::Integrator_BasisTimesScalar<EvaluationType,panzer::Traits>( plist ) );

	field_manager.registerEvaluator<EvaluationType>( source_op );
    }

    // Build the overall residual for the Poisson equation by summing the
    // above residuals.
    {
	Teuchos::RCP<std::vector<std::string> > residuals_to_sum = 
	    Teuchos::rcp( new std::vector<std::string> );

	if ( this->m_build_transient_support )
	{
	    residuals_to_sum->push_back( "RESIDUAL_TEMPERATURE_TRANSIENT_OP" );
	}
	residuals_to_sum->push_back( "RESIDUAL_TEMPERATURE_DIFFUSION_OP" );
	residuals_to_sum->push_back( "RESIDUAL_TEMPERATURE_SOURCE_OP" );

	Teuchos::ParameterList plist( "Temperature Residual" );
	plist.set( "Sum Name", "RESIDUAL_TEMPERATURE" );
	plist.set( "Values Names", residuals_to_sum );
	plist.set( "Data Layout", this->m_basis->functional );

	Teuchos::RCP< PHX::Evaluator<panzer::Traits> > residual_op =
	    Teuchos::rcp(
		new panzer::Sum<EvaluationType,panzer::Traits>( plist ) );

	field_manager.registerEvaluator<EvaluationType>( residual_op );
    }
}

} // end namespace TransientPoisson
} // end namespace Chimera

#endif // end TRANSIENTPOISSON_EQSET_IMPL_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_EquationSet_impl.hpp
//---------------------------------------------------------------------------//


