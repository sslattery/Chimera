//---------------------------------------------------------------------------//
// \file Chimera_OperatorTools.hpp
// \author Stuart R. Slattery
// \brief OperatorTools declaration.
//---------------------------------------------------------------------------//

#ifndef Chimera_OPERATORTOOLS_HPP
#define Chimera_OPERATORTOOLS_HPP

#include <vector>

#include <Teuchos_RCP.hpp>

#include <Epetra_Operator.h>

namespace Chimera
{
namespace Solvers
{

class OperatorTools
{
  public:

    // Constructor.
    OperatorTools();

    // Compute the spectral radius of an operator.
    double static spectralRadius( const Teuchos::RCP<Epetra_Operator>& matrix );

    // Compute the stiffness ratio of the operator.
    double static stiffnessRatio( const Teuchos::RCP<Epetra_Operator>& matrix );
};

} // end namespace Solvers
} // end namespace Chimera

#endif // Chimera_OPERATORTOOLS_HPP

//---------------------------------------------------------------------------//
// end Chimera_OperatorTools.hpp
//---------------------------------------------------------------------------//

