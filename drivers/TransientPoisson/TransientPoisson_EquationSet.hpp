//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_EquationSet.hpp
 * \author Stuart R. Slattery
 * \brief Equation set declaration for transient Poisson problem.
 */
//---------------------------------------------------------------------------//

#ifndef TRANSIENTPOISSON_EQSET_HPP
#define TRANSIENTPOISSON_EQSET_HPP

#include <vector>
#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Phalanx_FieldManager.hpp>

#include <Panzer_EquationSet_DefaultImpl.hpp>
#include <Panzer_Traits.hpp>
#include <Panzer_GlobalData.hpp>
#include <Panzer_InputEquationSet.hpp>
#include <Panzer_CellData.hpp>
#include <Panzer_BasisIRLayout.hpp>

//---------------------------------------------------------------------------//

namespace Chimera
{
namespace TransientPoisson
{

template<typename EvaluationType>
class TransientPoissonEquationSet 
    : public panzer::EquationSet_DefaultImpl<EvaluationType>
{
  public:

    // Constructor.
    TransientPoissonEquationSet( 
	const panzer::InputEquationSet& ies,
	const panzer::CellData& cell_data,
	const Teuchos::RCP<panzer::GlobalData>& global_data,
	const bool build_transient_support );

    // Build and register the evaluators with the field manager.
    void buildAndRegisterEquationSetEvaluators(
	PHX::FieldManager<panzer::Traits>& field_manager,
	const std::vector<std::pair<std::string,Teuchos::RCP<panzer::BasisIRLayout> > > & dofs,
	const Teuchos::ParameterList& user_data ) const;
};

} // end namespace TransientPoisson
} // end namespace Chimera

//---------------------------------------------------------------------------//
// Template includes.
//---------------------------------------------------------------------------//

#include "TransientPoisson_EquationSet_impl.hpp"

//---------------------------------------------------------------------------//

#endif // end TRANSIENTPOISSON_EQSET_HPP

//---------------------------------------------------------------------------//
// end TransientPoisson_EquationSet.hpp
//---------------------------------------------------------------------------//

