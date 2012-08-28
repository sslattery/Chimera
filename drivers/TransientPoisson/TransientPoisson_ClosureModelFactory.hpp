//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_ClosureModelFactory.hpp
 * \author Stuart R. Slattery
 * \brief Closure model factory for the transient Poisson problem.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_CLOSUREMODEL_HPP
#define TRANSIENTPOISSON_CLOSUREMODEL_HPP

#include <string>
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParamterList.hpp>

#include <Phalanx_Evaluator.hpp>
#include <Phalanx_FieldManager.hpp>

#include <Panzer_ClosureModel_Factory.hpp>
#include <Panzer_InputEquationSet.hpp>
#include <Panzer_Traits.hpp>

//---------------------------------------------------------------------------//

namespace Chimera
{
namespace TransientPoisson
{

template<typename EvaluationType>
class ClosureModelFactory : public panzer::ClosureModelFactory<EvaluationType>
{
  public:

    // Factory method.
   Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
   buildClosureModels( const std::string& model_id,
		       const panzer::InputEquationSet& equation_set,
		       const Teuchos::ParameterList& models,
		       const Teuchos::ParameterList& default_params,
		       const Teuchos::ParameterList& user_data,
		       const Teuchos::RCP<panzer::GlobalData>& global_data,
		       PHX::FieldManager<panzer::Traits>& field_manager ) const;

};

} // end namepsace TransientPoisson
} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "TransientPoisson_ClosureModelFactory_impl.hpp"

#endif // end TRANSIENTPOISSON_CLOSUREMODEL_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_ClosureModelFactory.hpp
//---------------------------------------------------------------------------//

