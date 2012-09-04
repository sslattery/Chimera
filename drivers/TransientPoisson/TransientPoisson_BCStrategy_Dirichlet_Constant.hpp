//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_BCStrategy_Dirichlet_Constant.hpp
 * \author Stuart R. Slattery
 * \brief Constant Dirichlet boundary condition declaration.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_BCSTRATEGY_DIRICHLET_CONSTANT_HPP
#define TRANSIENTPOISSON_BCSTRATEGY_DIRICHLET_CONSTANT_HPP

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Panzer_BCStrategy_Dirichlet_DefaultImpl.hpp>
#include <Panzer_Traits.hpp>
#include <Panzer_PureBasis.hpp>
#include <Panzer_PhysicsBlock.hpp>

#include "Phalanx_FieldManager.hpp"

namespace Chimera
{
namespace TransientPoisson
{

template<typename EvaluationType>
class BCStrategy_Dirichlet_Constant :
	public panzer::BCStrategy_Dirichlet_DefaultImpl<EvaluationType>
{
  public:
    
    BCStrategy_Dirichlet_Constant( const panzer::BC& bc,
				   const Teuchos::RCP<panzer::GlobalData>& global_data );

    void setup( const panzer::PhysicsBlock& side_pb,
		const Teuchos::ParameterList& user_data );

    void buildAndRegisterEvaluators( PHX::FieldManager<panzer::Traits>& fm,
				     cont panzer::PhysicsBlock& pb,
				     const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& factory,
				     const Teuchos::ParameterList& models,
				     const Teuchos::ParameterList& user_data ) const;

    // Data.
    std::string residual_name;
    Teuchos::RCP<panzer::PureBasis> basis;
};

} // end namespace TransientPoisson
} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "TransientPoisson_BCStrategy_Dirichlet_Constant_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end TRANSIENTPOISSON_BCSTRATEGY_DIRICHLET_CONSTANT_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_BCStrategy_Dirichlet_Constant.hpp
//---------------------------------------------------------------------------//

