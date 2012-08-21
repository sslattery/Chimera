//---------------------------------------------------------------------------//
/*!
 * \file IterativeSolver.hpp
 * \author Stuart R. Slattery
 * \brief Protocol declaration for iterative solvers.
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_ITERATIVESOLVER_HPP
#define CHIMERA_ITERATIVESOLVER_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_ParameterList.hpp>

#include <BelosLinearProblem.hpp>

namespace Chimera
{

/*!
 * \brief Protocol definition for iterative solvers.
 *
 * Scalar must implement Teuchos::ScalarTraits, MV must implement
 * Belos::MultiVecTraits, and OP must implement Belos::OperatorTraits.
 */
template<typename Scalar, typename MV, typename OP>
class IterativeSolver : public Teuchos::Describable
{
  public:

    //@{
    //! Typedefs.
    typedef Scalar                                 scalar_type;
    typedef MV                                     multivector_type;
    typedef OP                                     operator_type;
    typedef Belos::LinearProblem<Scalar,MV,OP>     LinearProblem;
    //@}

    /*!
     * \brief Default constructor.
     */
    IterativeSolver()
    { /* ... */ }

    /*!
     * \brief Destructor.
     */
    virtual ~IterativeSolver()
    { /* ... */ }

    /*!
     * \brief Iterate to convergence.
     */
    virtual void iterate() = 0;

    /*!
     * \brief Get the iteration count from the last solve.
     */
    virtual std::size_t getNumIters() = 0;

    /*!
     * \brief Get the convergence tolerance.
     */
    virtual double getTolerance() = 0;

    /*!
     * \brief Get the iteration number limit.
     */
    virtual std::size_t getMaxIters() = 0;

    /*!
     * \brief Get the valid parameters for the solver.
     */
    virtual const Teuchos::ParamterList& getValidParameters() = 0;

    /*!
     * \brief Get the current solver parameters.
     */
    virtual const Teuchos::ParameterList& getParameterList() = 0;

    /*! 
     * \brief Get the linear problem.
     */
    virtual const LinearProblem& getLinearProblem() = 0;
};

} // end namespace Chimera

#end // CHIMERA_ITERATIVESOLVER_HPP

//---------------------------------------------------------------------------//
// end IterativeSolver.hpp
//---------------------------------------------------------------------------//

