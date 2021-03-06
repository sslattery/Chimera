//---------------------------------------------------------------------------//
/*! 
 * \file VTKOutput.hpp
 * \author Stuart R. Slattery
 * \brief VTK file output declaration for mesh and fields.
 */
//---------------------------------------------------------------------------//

#ifndef CHIMERA_VTKOUTPUT_HPP
#define CHIMERA_VTKOUTPUT_HPP

#include <string>
#include <fstream>

#include <Partitioner.hpp>

#include <Teuchos_RCP.hpp>

#include <Tpetra_Vector.hpp>

namespace Chimera
{

class VTKOutput
{

  public:
    
    //@{
    //! Typedefs.
    typedef Partitioner::VecPair             VecPair;
    typedef Partitioner::RCP_Comm            RCP_Comm;
    typedef Partitioner::RCP_ParameterList   RCP_ParameterList;
    typedef Teuchos::RCP<Partitioner>        RCP_Partitioner;
    typedef Tpetra::Vector<double>           VectorType;
    typedef Teuchos::RCP<VectorType>         RCP_Vector;
    //@}

    //! Field topology enum.
    enum FieldTopology{ VERTEX_FIELD = 0, CELL_FIELD };

    // Constructor.
    VTKOutput( const RCP_Comm &comm, 
	       const RCP_Partitioner &partitioner,
	       const RCP_ParameterList &plist );

    // Destructor.
    ~VTKOutput();

    // Add a field to the file.
    void addField( const int field_type,
		   const RCP_Vector &field,
		   const std::string &name );

    // Write the file.
    void write();

  private:

    // Communicator.
    RCP_Comm d_comm;

    // Partitioner.
    RCP_Partitioner d_partitioner;

    // Mesh file.
    std::ofstream d_vtk_file;

    // Number of vertices in the x_direction.
    std::size_t d_nx;
    
    // Number of vertices in the y_direction.
    std::size_t d_ny;
};

} // end namespace Chimera

#endif // CHIMERA_VTKOUTPUT_HPP

//---------------------------------------------------------------------------//
// end VTKOutput.hpp
//---------------------------------------------------------------------------//

