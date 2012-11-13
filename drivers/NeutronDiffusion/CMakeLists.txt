INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Diffusion
  DEPLIBS chimera
  SOURCES Partitioner.cpp DiffusionProblem.cpp Mesh.cpp VTKOutput.cpp neutron_diffusion.cpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  )
