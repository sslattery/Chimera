//---------------------------------------------------------------------------//
/*!
 * \file DiffusionProblem.hpp
 * \author Stuart R. Slattery
 * \brief Diffusion problem declaration.
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_DIFFUSIONPROBLEM_HPP
#define CHIMERA_DIFFUSIONPROBLEM_HPP

#include "Partitioner.hpp"

#include <Chimera_LinearProblem.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Chimera
{

class DiffusionProblem
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<Partitioner>                   RCP_Partitioner;
    typedef LinearProblem<double,int,int>               LinearProblemType;
    typedef Teuchos::RCP<LinearProblemType>             RCP_LinearProblem;
    typedef Teuchos::RCP<Teuchos::ParameterList>        RCP_ParameterList;
    typedef Teuchos::Comm<int>                          CommType;
    typedef Teuchos::RCP<const CommType>                RCP_Comm;
    //@}

    // Constructor.
    DiffusionProblem( const RCP_Comm& comm, 
		      const RCP_Partitioner& partitioner,
		      const RCP_ParameterList& plist,
		      bool jacobi_precondition = false );

    // Destructor.
    ~DiffusionProblem();

    // Get the diffusion problem.
    RCP_LinearProblem getProblem() 
    { return d_linear_problem; }

  private:

    // The diffusion problem.
    RCP_LinearProblem d_linear_problem;
};

} // end namespace Chimera

#endif // end CHIMERA_DIFFUSIONPROBLEM_HPP

//---------------------------------------------------------------------------//
// end DiffusionProblem.hpp
//---------------------------------------------------------------------------//

