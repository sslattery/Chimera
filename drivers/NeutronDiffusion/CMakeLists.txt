INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Diffusion
  DEPLIBS belos chimera mcls
  SOURCES Partitioner.cpp DiffusionProblem.cpp Mesh.cpp VTKOutput.cpp neutron_diffusion.cpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(
  DiffusionINP
  SOURCE_FILES input.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS Diffusion
)
