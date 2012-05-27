//---------------------------------------------------------------------------//
/*! 
 * \file VTKWriter.hpp
 * \author Stuart R. Slattery
 * \brief VTK file writer declaration for mesh and fields.
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_VTKWRITER_HPP
#define CHIMERA_VTKWRITER_HPP

#include <fstream>

#include <Partitioner.hpp>

#include <Teuchos_RCP.hpp>

#include <Tpetra_Vector.hpp>

namespace Chimera
{

class VTKWriter
{

  public:
    
    //@{
    //! Typedefs.
    typedef Partitioner::VecPair             VecPair;
    typedef Partitioner::RCP_Comm            RCP_Comm;
    typedef Partitioner::RCP_ParameterList   RCP_ParameterList;
    typedef Teuchos::RCP<Partitioner>        RCP_Partitioner;
    //@}

    // Constructor.
    VTKWriter( const RCP_Comm &comm, 
	       const RCP_Partitioner &partitioner,
	       const RCP_ParameterList &plist );

    // Destructor.
    ~VTKWriter();

    // Add a field to the file.
    void addField( const Tpetra::Vector<double>& field );

    // Write the file.
    void write();

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Mesh file.
    std::ofstream d_vtk_file;
};

} // end namespace Chimera

#endif // CHIMERA_VTKWRITER_HPP

//---------------------------------------------------------------------------//
// end VTKWriter.hpp
//---------------------------------------------------------------------------//

