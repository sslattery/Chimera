//---------------------------------------------------------------------------//
/*!
 * \file VTKOutput.cpp
 * \author Stuart R. Slattery
 * \brief VTK file output definition for mesh and fields.
 */
//---------------------------------------------------------------------------//

#include <vector>
#include <sstream>

#include "VTKOutput.hpp"

#include <Chimera_Assertion.hpp>

#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_CommHelpers.hpp>

namespace Chimera
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
VTKOutput::VTKOutput( const RCP_Comm &comm, 
		      const RCP_Partitioner &partitioner,
		      const RCP_ParameterList &plist )
    : d_comm( comm )
    , d_mesh( partitioner->getMesh() )
{
    // Local mesh size.
    d_nx = partitioner->getMesh()->getLocalEdges().first.size();
    d_ny = partitioner->getMesh()->getLocalEdges().second.size();

    // Open the file to write.
    std::stringstream convert_rank;
    convert_rank << d_comm->getRank();
    std::stringstream convert_size;
    convert_size << d_comm->getSize();

    std::string filename = 
	plist->get<std::string>( "OUTPUT_FILENAME" ) 
	+ "." + convert_size.str() + "." + convert_rank.str();

    d_vtk_file.open( filename.c_str() );

    // File header.
    d_vtk_file << "# vtk DataFile Version 3.0\n";
    d_vtk_file << "Chimera version 0.1\n";
    d_vtk_file << "ASCII\n";
    d_vtk_file << "DATASET RECTILINEAR_GRID\n";
    d_vtk_file << "DIMENSIONS " << d_nx << " " << d_ny << " 1\n";

    // Set the mesh.
    std::vector<double>::const_iterator vec_it;

    d_vtk_file << "X_COORDINATES " << d_nx << " double\n";
    for ( vec_it = partitioner->getMesh()->getLocalEdges().first.begin();
	  vec_it != partitioner->getMesh()->getLocalEdges().first.end();
	  ++vec_it )
    {
	d_vtk_file << *vec_it << " ";
    }
    d_vtk_file << "\n";

    d_vtk_file << "Y_COORDINATES " << d_ny << " double\n";
    for ( vec_it = partitioner->getMesh()->getLocalEdges().second.begin();
	  vec_it != partitioner->getMesh()->getLocalEdges().second.end();
	  ++vec_it )
    {
	d_vtk_file << *vec_it << " ";
    }
    d_vtk_file << "\n";

    d_vtk_file << "Z_COORDINATES 1 double\n";
    d_vtk_file << "0\n";
}

//---------------------------------------------------------------------------//
/*!
 * \brief Destructor.
 */
VTKOutput::~VTKOutput()
{ /* ... */ }

//---------------------------------------------------------------------------//
/*!
 * \brief Add a field to the file.
 */
void VTKOutput::addField( const int field_type,
			  const RCP_Vector &field,
			  const std::string &name )
{
    // Setup.
    int num_components = 1;
    int local_length = field->getLocalLength();
    Teuchos::ArrayView<const int> local_rows = 
	field->getMap()->getNodeElementList();

    // Field header.
    if ( field_type == VERTEX_FIELD )
    {
	d_vtk_file << "POINT_DATA " << local_length << "\n";
    }
    else if ( field_type == CELL_FIELD )
    {
	d_vtk_file << "CELL_DATA " << local_length << "\n";
    }

    d_vtk_file << "SCALARS " << name << " double " << num_components << "\n";
    d_vtk_file << "LOOKUP_TABLE default\n";

    // Write to file.
    for ( int i = 0; i < local_length; ++i )
    {
	d_vtk_file << field->get1dView()[i] << "\n";
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write the file.
 */
void VTKOutput::write()
{
    d_vtk_file.close();
}

//---------------------------------------------------------------------------//

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end VTKOutput.cpp
//---------------------------------------------------------------------------//

