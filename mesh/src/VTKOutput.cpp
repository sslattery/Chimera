//---------------------------------------------------------------------------//
/*!
 * \file VTKOutput.cpp
 * \author Stuart R. Slattery
 * \brief VTK file output definition for mesh and fields.
 */
//---------------------------------------------------------------------------//

#include <vector>

#include "VTKOutput.hpp"

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
{
    if ( d_comm->getRank() == 0 )
    {
	// Open the file to write.
	d_vtk_file.open( plist->get<std::string>( "OUTPUT_FILENAME" ).c_str() );

	// File header.
	d_vtk_file << "# vtk DataFile Version 2.0\n";
	d_vtk_file << "Chimera version 0.1\n";
	d_vtk_file << "ASCII\n";
	d_vtk_file << "DATASET RECTILINEAR_GRID\n";

	// Set the mesh.
	std::vector<double>::const_iterator vec_it;
	d_nx = partitioner->getGlobalEdges().first.size();
	d_ny = partitioner->getGlobalEdges().second.size();
	int nz = 1;
	d_vtk_file << "DIMENSIONS " << d_nx << " " << d_ny << " " << nz << "\n";

	d_vtk_file << "X_COORDINATES " << d_nx << " double\n";
	for ( vec_it = partitioner->getGlobalEdges().first.begin();
	      vec_it != partitioner->getGlobalEdges().first.end();
	      ++vec_it )
	{
	    d_vtk_file << *vec_it << " ";
	}
	d_vtk_file << "\n";

	d_vtk_file << "Y_COORDINATES " << d_ny << " double\n";
	for ( vec_it = partitioner->getGlobalEdges().second.begin();
	      vec_it != partitioner->getGlobalEdges().second.end();
	      ++vec_it )
	{
	    d_vtk_file << *vec_it << " ";
	}
	d_vtk_file << "\n";

	d_vtk_file << "Z_COORDINATES " << nz << " double\n";
	d_vtk_file << "0\n";
	
	d_vtk_file << "0\n"
    }
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
void VTKOutput::addField( const int field_topology,
			  const int field_type,
			  const RCP_MultiVector &field, 
			  const std::string &name )
{
    // Check the input arguments.
    int global_num_verts = d_nx*d_ny;
    int global_num_cells = (d_nx-1)*(d_ny-1);
    int num_components = field->getNumVectors();
    int global_length = field->getGlobalLength();

    testPrecondition( num_components < 3, "Number of field components > 2" );

    if ( field_topology == VERTEX_FIELD )
    {
	testPrecondition( global_length == global_num_verts,
			  "Field size != number of vertices" );
    }
    else if ( field_topology == CELL_FIELD )
    {
	testPrecondition( global_length == global_num_cells,
			  "Field size != number of cells" );
    }
    else
    {
	testPrecondition( field_topology == VERTEX_FIELD || 
			  field_topology == CELL_FIELD,
			  "Field topology not supported" );
    }

    testPrecondition( field_type == SCALAR_FIELD ||
		      field_type == VECTOR_FIELD,
		      "Field type not supported" );

    // Field header.
    int my_rank = d_comm->getRank();
    if ( my_rank == 0 )
    {
	if ( field_topology == VERTEX_FIELD )
	{
	    d_vtk_file << "POINT_DATA " << global_num_verts << "\n";
	}
	else if ( field_topology == CELL_FIELD )
	{
	    d_vtk_file << "CELL_DATA " << global_num_cells << "\n";
	}

	if ( field_type == SCALAR_FIELD )
	{
	    d_vtk_file << "SCALARS " << name << " double" 
		       << num_components << "\n";
	}
	else if ( field_type == VECTOR_FIELD )
	{
	    d_vtk_file << "VECTORS " << name << " double" << "\n";
	}

	d_vtk_file << "LOOKUP TABLE default\n";
    }

    // Serialize each field component.
    for ( int n = 0; n < num_components )
    {

    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Write the file.
 */
void VTKOutput::write()
{
    if ( d_comm->getRank() == 0 )
    {
	d_vtk_file.close();
    }
}

} // end namespace Chimera

//---------------------------------------------------------------------------//
// end VTKOutput.cpp
//---------------------------------------------------------------------------//

