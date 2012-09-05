//---------------------------------------------------------------------------//
// \file OperatorTools.hpp
// \author Stuart R. Slattery
// \brief OperatorTools declaration.
//---------------------------------------------------------------------------//

#ifndef HMCSA_OPERATORTOOLS_HPP
#define HMCSA_OPERATORTOOLS_HPP

#include <vector>

#include <Teuchos_RCP.hpp>

#include <Epetra_Operator.h>

namespace HMCSA
{

namespace OperatorTools
{

// Compute the spectral radius of an operator.
double spectralRadius( const Teuchos::RCP<Epetra_Operator>& matrix );

// Compute the stiffness ratio of the operator.
double stiffnessRatio( const Teuchos::RCP<Epetra_Operator>& matrix );

} // end namespace OperatorTools

} // end namespace HMCSA

#endif // HMCSA_OPERATORTOOLS_HPP

//---------------------------------------------------------------------------//
// end OperatorTools.hpp
//---------------------------------------------------------------------------//

