//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_ClosureModelFactory_impl.hpp
 * \author Stuart R. Slattery
 * \brief Closure model factory for the transient Poisson problem.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_CLOSUREMODEL_IMPL_HPP
#define TRANSIENTPOISSON_CLOSUREMODEL_IMPL_HPP

#include <Panzer_Constant.hpp>

namespace Chimera
{
namespace TransientPoisson
{
//---------------------------------------------------------------------------//
/*!
 * \brief Closure model factory method.
 */
template<typename EvaluationType>
Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
ClosureModelFactory<EvaluationType>::buildClosureModels( 
    const std::string& model_id,
    const panzer::InputEquationSet& equation_set,
    const Teuchos::ParameterList& models,
    const Teuchos::ParameterList& default_params,
    const Teuchos::ParameterList& user_data,
    const Teuchos::RCP<panzer::GlobalData>& global_data,
    PHX::FieldManager<panzer::Traits>& field_manager ) const
{
    // Setup a vector for the closure evaluators.
    Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >  =
	Teuchos::rcp( 
	    new std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > );

    // Check that this model id is valid.
    if ( !models.isSublist( model_id ) )
    {
	models.print( std::cout );
	std::stringstream msg;
	msg << "Falied to find requested model, \"" << model_id 
	    << "\", for equation set:\n" << std::endl;
	TEUCHOS_TEST_FOR_EXCEPTION( !models.isSublist(model_id), 
				    std::logic_error, msg.str() );
    }

    // Loop through the models and build evaluators.
    const Teuchos::ParameterList& my_models = models.sublist( model_id );

    for ( Teuchos::ParameterList::ConstIterator model_it = my_models.begin();
	  model_it != my_models.end(); ++model_it )
    {
	bool found = false;

	const std::string key = model_it->first;
	Teuchos::ParameterList input;
	const Teuchos::ParameterEntry& entry = model_it->second;
	const Teuchos::ParameterList& plist = 
	    Teuchos::getValue<Teuchos::ParameterList>( entry );

	// Create evaluators of type double.
	if ( plist.isType<double>( "Value" ) )
	{
	    // Set constant scalar evaluator at integration points.
	    {
		input.set( "Name", key );
		input.set( "Value", plist.get<double>( "Value" ) );
		input.set( "Data Layout", 
			   default_params.get<Teuchos::RCP<panzer::IntegrationRule> >( "IR" )->dl_scalar );
	      Teuchos:RCP<PHX::Evaluator<EvaluationType,panzer::Traits> > constant_evaluator =
		    Teuchos::rcp( new panzer::Constant<EvaluationType,panzer::Traits>( input ) );
		evaluators.push_back( constant_evaluator );
	    }

	    // Set constant scalar evaluator at basis supports.
	    {
		input.set( "Name", key );
		input.set( "Value", plist.get<double>( "Value" ) );
		input.set( "Data Layout", 
			   default_params.get<Teuchos::RCP<panzer::BasisIRLayout> >( "Basis" )->functional );
	      Teuchos:RCP<PHX::Evaluator<EvaluationType,panzer::Traits> > constant_evaluator =
		    Teuchos::rcp( new panzer::Constant<EvaluationType,panzer::Traits>( input ) );
		evaluators.push_back( constant_evaluator );
	    }

	    found = true;
	}

	// Check that the closure for this model was found.
	if (!found) {
	    p      std::stringstream msg;
	    msg << "ClosureModelFactory failed to build evaluator for key \"" << key 
		<< "\"\nin model \"" << model_id 
		<< "\".  Please correct the type or add support to the \nfactory." <<std::endl;
	    TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, msg.str());
	}
    }

    return evaluators;
}

//---------------------------------------------------------------------------//

} // end namespace TransientPoisson
} // end namespace Chimera

#endif // end TRANSIENTPOISSON_CLOSUREMODEL_IMPL_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_ClosureModelFactory_impl.hpp
//---------------------------------------------------------------------------//

