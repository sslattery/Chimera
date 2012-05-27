//---------------------------------------------------------------------------//
/*!
 * \file Partioner.hpp
 * \author Stuart R. Slattery
 * \brief Mesh partitioner declaration.
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_PARTITIONER_HPP
#define CHIMERA_PARTITIONER_HPP

#include "Mesh.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

namespace Chimera
{

class Partitioner
{
  public:

    //@{
    //! typedefs.
    typedef Mesh::SizePair                        SizePair;
    typedef Teuchos::RCP<Mesh>                    RCP_Mesh;
    typedef Teuchos::Comm<int>                    CommType;
    typedef Teuchos::RCP<const CommType>          RCP_Comm;
    typedef Teuchos::RCP<Teuchos::ParameterList>  RCP_ParameterList;
    //@}

    // Constructor.
    Partitioner( const RCP_Comm &comm, const RCP_ParameterList &plist );

    // Destructor.
    ~Partitioner();

    // Get the number of blocks.
    const SizePair& getNumBlocks() const
    { return d_num_blocks; }

    // Get the mesh.
    const RCP_Mesh& getMesh() const
    { return d_mesh; }

  private:

    // Number of blocks.
    SizePair d_num_blocks;

    // Mesh.
    RCP_Mesh d_mesh;
};

} // end namespace Chimera

#endif // end CHIMERA_PARTITIONER_HPP

//---------------------------------------------------------------------------//
// end Partitioner.hpp
//---------------------------------------------------------------------------//

