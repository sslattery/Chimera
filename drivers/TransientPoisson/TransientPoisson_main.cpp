//---------------------------------------------------------------------------//
/*!
 * \file TransientPoisson_main.cpp
 * \author Stuart R. Slattery
 * \brief Driver for transient poisson problem.
 */
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <string>
#include <set>

#include "Chimera_MCSA.hpp"
#include "Chimera_JacobiPreconditioner.hpp"

#include "TransientPoisson_EquationSetFactory.hpp"
#include "TransientPoisson_ClosureModelFactory_TemplateBuilder.hpp"
#include "TransientPoisson_BCStrategy_Factory.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_MpiComm.h>
#include <Epetra_LinearProblem.h>

#include <Panzer_config.hpp>
#include <Panzer_InputPhysicsBlock.hpp>
#include <Panzer_BC.hpp>
#include <Panzer_PhysicsBlock.hpp>
#include <Panzer_InputEquationSet.hpp>
#include <Panzer_CellData.hpp>
#include <Panzer_PureBasis.hpp>
#include <Panzer_Workset_Builder.hpp>
#include <Panzer_WorksetContainer.hpp>
#include <Panzer_ParameterList_ObjectBuilders.hpp>
#include <Panzer_AssemblyEngine.hpp>
#include <Panzer_AssemblyEngine_TemplateManager.hpp>
#include <Panzer_AssemblyEngine_TemplateBuilder.hpp>
#include <Panzer_LinearObjFactory.hpp>
#include <Panzer_EpetraLinearObjFactory.hpp>
#include <Panzer_DOFManagerFactory.hpp>
#include <Panzer_DOFManager.hpp>
#include <Panzer_FieldManagerBuilder.hpp>

#include <Panzer_STK_config.hpp>
#include <Panzer_STK_Interface.hpp>
#include <Panzer_STK_SquareQuadMeshFactory.hpp>
#include <Panzer_STK_WorksetFactory.hpp>
#include <Panzer_STKConnManager.hpp>
#include <Panzer_STK_SetupUtilities.hpp>
#include <Panzer_STK_Utilities.hpp>
#include <Panzer_STK_Version.hpp>

#include <AztecOO.h>

int main( int argc, char * argv[] )
{
    // Initialize parallel communication.
    Teuchos::GlobalMPISession mpi_session( &argc, &argv );
    Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp( 
	new Epetra_MpiComm( MPI_COMM_WORLD ) );
    Teuchos::FancyOStream out( Teuchos::rcpFromRef( std::cout ) );
    out.setOutputToRootOnly( 0 );
    out.setShowProcRank( true );

    // Create transient poisson factories.
    Chimera::TransientPoisson::EquationSetFactory eqset_factory;
    Chimera::TransientPoisson::BCStrategyFactory bc_factory;
    Chimera::TransientPoisson::ClosureModelFactory_TemplateBuilder cm_builder;

    // Generate the mesh.
    const std::size_t workset_size = 20;
    panzer_stk::SquareQuadMeshFactory mesh_factory;

    Teuchos::RCP<Teuchos::ParameterList> plist = 
	Teuchos::rcp( new Teuchos::ParameterList );
    plist->set( "X Blocks", 1 );
    plist->set( "Y Blocks", 1 );
    plist->set( "X Elements", 100 );
    plist->set( "Y Elements", 100 );
    mesh_factory.setParameterList( plist );

    Teuchos::RCP<panzer_stk::STK_Interface> mesh = 
	mesh_factory.buildUncommitedMesh( MPI_COMM_WORLD );

    // Build the physics blocks.
    panzer::InputPhysicsBlock input_blocks;
    std::vector<panzer::BC> boundary_conditions;
    std::vector<Teuchos::RCP<panzer::PhysicsBlock> > physics_blocks;
    {
	// Switch for transient support.
	bool build_transient_support = false;
	
	// Build the input equation set and block.
	panzer::InputEquationSet input_equations;
	input_equations.name = "TransientPoisson";
	input_equations.basis = "Q1";
	input_equations.integration_order = 2;
	input_equations.model_id = "solid";
	input_equations.prefix = "";
	input_blocks.physics_block_id = "4";
	input_blocks.eq_sets.push_back( input_equations );

	// Left side boundary condition.
	{
	    std::size_t bc_id = 0;
	    panzer::BCType bc_type = panzer::BCT_Dirichlet;
	    std::string sideset_id = "left";
	    std::string element_block_id = "eblock-0_0";
	    std::string dof_name = "TEMPERATURE";
	    std::string strategy = "Constant";
	    double value = 5.0;
	    Teuchos::ParameterList bc_list;
	    bc_list.set( "Value", value );
	    panzer::BC bc( bc_id, bc_type, sideset_id, element_block_id,
			   dof_name, strategy, bc_list );
	    boundary_conditions.push_back( bc );
	}

	// Top side boundary condition.
	{
	    std::size_t bc_id = 1;
	    panzer::BCType bc_type = panzer::BCT_Dirichlet;
	    std::string sideset_id = "top";
	    std::string element_block_id = "eblock-0_0";
	    std::string dof_name = "TEMPERATURE";
	    std::string strategy = "Constant";
	    double value = 0.0;
	    Teuchos::ParameterList bc_list;
	    bc_list.set( "Value", value );
	    panzer::BC bc( bc_id, bc_type, sideset_id, element_block_id,
			   dof_name, strategy, bc_list );
	    boundary_conditions.push_back( bc );
	}

	// Right side boundary condition.
	{
	    std::size_t bc_id = 2;
	    panzer::BCType bc_type = panzer::BCT_Dirichlet;
	    std::string sideset_id = "right";
	    std::string element_block_id = "eblock-0_0";
	    std::string dof_name = "TEMPERATURE";
	    std::string strategy = "Constant";
	    double value = 5.0;
	    Teuchos::ParameterList bc_list;
	    bc_list.set( "Value", value );
	    panzer::BC bc( bc_id, bc_type, sideset_id, element_block_id,
			   dof_name, strategy, bc_list );
	    boundary_conditions.push_back( bc );
	}

	// Bottom side boundary condition.
	{
	    std::size_t bc_id = 3;
	    panzer::BCType bc_type = panzer::BCT_Dirichlet;
	    std::string sideset_id = "bottom";
	    std::string element_block_id = "eblock-0_0";
	    std::string dof_name = "TEMPERATURE";
	    std::string strategy = "Constant";
	    double value = 0.0;
	    Teuchos::ParameterList bc_list;
	    bc_list.set( "Value", value );
	    panzer::BC bc( bc_id, bc_type, sideset_id, element_block_id,
			   dof_name, strategy, bc_list );
	    boundary_conditions.push_back( bc );
	}

	// Get the cell data.
	const panzer::CellData volume_cell_data(
	    workset_size,
	    mesh->getCellTopology( "eblock-0_0" )->getDimension(),
	    mesh->getCellTopology( "eblock-0_0" ) );

	// Build global data.
	Teuchos::RCP<panzer::GlobalData> global_data = 
	    panzer::createGlobalData();

	// Generate the physics block.
	Teuchos::RCP<panzer::PhysicsBlock> block = Teuchos::rcp(
	    new panzer::PhysicsBlock( input_blocks, "eblock-0_0",
				      volume_cell_data, eqset_factory,
				      global_data, build_transient_support ) );
	physics_blocks.push_back( block );
    }

    // Finalize the mesh by adding DOFs.
    {
	Teuchos::RCP<panzer::PhysicsBlock> block = physics_blocks[0];
	const std::vector<panzer::StrPureBasisPair>& block_fields =
	    block->getProvidedDOFs();

	std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp> field_names;
	field_names.insert( block_fields.begin(), block_fields.end() );

	std::set<panzer::StrPureBasisPair,panzer::StrPureBasisComp>::const_iterator
	    field_name_iterator;
	for ( field_name_iterator = field_names.begin();
	      field_name_iterator != field_names.end();
	      ++field_name_iterator )
	{
	    mesh->addSolutionField( field_name_iterator->first,
				    block->elementBlockID() );
	}

	mesh_factory.completeMeshConstruction( *mesh, MPI_COMM_WORLD );
    }

    // Build the worksets.
    Teuchos::RCP<panzer_stk::WorksetFactory> workset_factory = Teuchos::rcp(
	new panzer_stk::WorksetFactory( mesh ) );
    Teuchos::RCP<panzer::WorksetContainer> workset_container = Teuchos::rcp(
	new panzer::WorksetContainer( workset_factory, physics_blocks,
				      workset_size ) );

    // Build the DOF manager.
    const Teuchos::RCP<panzer::ConnManager<int,int> > conn_manager =
	Teuchos::rcp( new panzer_stk::STKConnManager( mesh ) );

    panzer::DOFManagerFactory<int,int> global_indexer_factory;
    Teuchos::RCP<panzer::UniqueGlobalIndexer<int,int> > dof_manager =
	global_indexer_factory.buildUniqueGlobalIndexer( 
	    Teuchos::opaqueWrapper(MPI_COMM_WORLD),
	    physics_blocks,
	    conn_manager );

    // Build the linear algebra object factory.
    Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > lin_obj_factory =
	Teuchos::rcp( new panzer::EpetraLinearObjFactory<panzer::Traits,int>(
			  comm.getConst(), dof_manager ) );

    // Setup the closure model for the source.
    panzer::ClosureModelFactory_TemplateManager<panzer::Traits> 
	closure_model_factory;
    closure_model_factory.buildObjects( cm_builder );

    Teuchos::ParameterList closure_models( "Closure Models" );
    closure_models.sublist( "solid" ).sublist( "SOURCE_TEMPERATURE" ).set<double>("Value", 0.0 );

    // Set user data.
    Teuchos::ParameterList user_data( "User Data" );
    user_data.set<double>("Thermal Conductivity", 2.0 );

    // Setup the field managers.
    Teuchos::RCP<panzer::FieldManagerBuilder<int,int> > field_manager_builder =
	Teuchos::rcp( new panzer::FieldManagerBuilder<int,int> );
    field_manager_builder->setupVolumeFieldManagers( *workset_container,
						     physics_blocks,
						     closure_model_factory,
						     closure_models,
						     *lin_obj_factory,
						     user_data );
    field_manager_builder->setupBCFieldManagers( *workset_container,
						 boundary_conditions,
						 physics_blocks,
						 eqset_factory,
						 closure_model_factory,
						 bc_factory,
						 closure_models,
						 *lin_obj_factory,
						 user_data );

    // Generate the assembly engine.
    panzer::AssemblyEngine_TemplateManager<panzer::Traits,int,int> assembly_engine;
    panzer::AssemblyEngine_TemplateBuilder<int,int> assembly_builder( 
	field_manager_builder, lin_obj_factory );
    assembly_engine.buildObjects( assembly_builder );

    // Assemble the linear system.
    Teuchos::RCP<panzer::LinearObjContainer> ghost_container =
	lin_obj_factory->buildGhostedLinearObjContainer();
    Teuchos::RCP<panzer::LinearObjContainer> container =
	lin_obj_factory->buildLinearObjContainer();
    lin_obj_factory->initializeGhostedContainer( panzer::LinearObjContainer::X |
						 panzer::LinearObjContainer::F |
						 panzer::LinearObjContainer::Mat |
						 panzer::LinearObjContainer::DxDt,
						 *ghost_container );
    lin_obj_factory->initializeContainer( panzer::LinearObjContainer::X |
					  panzer::LinearObjContainer::F |
					  panzer::LinearObjContainer::Mat |
					  panzer::LinearObjContainer::DxDt,
					  *container );
    ghost_container->initialize();
    container->initialize();

    // Evaluate the physics to build the residual and Jacobian.
    panzer::AssemblyEngineInArgs input( ghost_container, container );
    input.alpha = 0;
    input.beta = 1;
    assembly_engine.getAsObject<panzer::Traits::Jacobian>()->evaluate( input );

    // Solve the linear problem with MCSA.
    Teuchos::RCP<panzer::EpetraLinearObjContainer> ep_container =
	Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>( container );

    Teuchos::RCP<Epetra_LinearProblem> problem = Teuchos::rcp(
	new Epetra_LinearProblem( &*ep_container->get_A(),
				  &*ep_container->get_x(),
				  &*ep_container->get_f() ) );

    Chimera::Solvers::JacobiPreconditioner preconditioner( problem );
    preconditioner.precondition();
    
    Chimera::Solvers::MCSA solver( problem );
    int max_iters = 1000;
    double tolerance = 1.0e-8;
    int num_histories = 50;
    double weight_cutoff = 1.0e-4;
    solver.iterate( max_iters, tolerance, num_histories, weight_cutoff );

    ep_container->get_x()->Scale( -1.0 );

    // Generate output.
    lin_obj_factory->globalToGhostContainer( 
	*container, *ghost_container,
	panzer::EpetraLinearObjContainer::X | 
	panzer::EpetraLinearObjContainer::DxDt );

    Teuchos::RCP<panzer::EpetraLinearObjContainer> ep_ghost_container =
	Teuchos::rcp_dynamic_cast<panzer::EpetraLinearObjContainer>( ghost_container );
    panzer_stk::write_solution_data( 
	*Teuchos::rcp_dynamic_cast<panzer::DOFManager<int,int> >( dof_manager ),
	*mesh, *ep_ghost_container->get_x() );
    mesh->writeToExodus( "transient_poisson.exo" );

    // Complete.
    return 0;
}

//---------------------------------------------------------------------------//
// end TransientPoisson_main.cpp
//---------------------------------------------------------------------------//

