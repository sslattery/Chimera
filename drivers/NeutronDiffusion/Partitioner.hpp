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
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>

namespace Chimera
{

class Partitioner
{
  public:

    //@{
    //! Typedefs.
    typedef Mesh::SizePair                        SizePair;
    typedef Mesh::VecPair                         VecPair;
    typedef Teuchos::RCP<Mesh>                    RCP_Mesh;
    typedef Teuchos::Comm<int>                    CommType;
    typedef Teuchos::RCP<const CommType>          RCP_Comm;
    typedef Teuchos::RCP<Teuchos::ParameterList>  RCP_ParameterList;
    //@}

    // Constructor.
    Partitioner( const RCP_Comm &comm, const RCP_ParameterList &plist );

    // Destructor.
    ~Partitioner();

    //! Get the number of blocks.
    const SizePair& getNumBlocks() const
    { return d_num_blocks; }

    //! Get the global edge vectors.
    const VecPair& getGlobalEdges() const
    { return d_global_edges; }

    //! Get the mesh.
    const RCP_Mesh& getMesh() const
    { return d_mesh; }

    //! Get the local rows (vertex-based global ids).
    Teuchos::ArrayView<int> getLocalRows()
    { return d_local_rows(); }

    //! Get the cell sizes.
    std::pair<double,double> getCellSizes()
    { return d_cell_size; }

  private:

    // Number of blocks.
    SizePair d_num_blocks;

    // Global edge vectors.
    VecPair d_global_edges;

    // Mesh.
    RCP_Mesh d_mesh;

    // Cell size.
    std::pair<double,double> d_cell_size;

    // Local rows (vertex-based global ids).
    Teuchos::Array<int> d_local_rows;
};

} // end namespace Chimera

#endif // end CHIMERA_PARTITIONER_HPP

//---------------------------------------------------------------------------//
// end Partitioner.hpp
//---------------------------------------------------------------------------//

