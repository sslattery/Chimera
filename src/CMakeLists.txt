INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Chimera_AdjointMC.hpp
  Chimera_BoostRNG.hpp
  Chimera_JacobiPreconditioner.hpp
  Chimrea_LinearProblem.hpp
  Chimera_MCSA.hpp
  Chimera_OperatorTools.hpp
  Chimera_RNGTraits.hpp
  ) 

APPEND_SET(SOURCES
  Chimera_AdjointMC.cpp
  Chimera_JacobiPreconditioner.cpp
  Chimera_MCSA.cpp
  Chimera_OperatorTools.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  chimera
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

INSTALL(FILES 
  ${HEADERS} 
  DESTINATION 
  ${CMAKE_INSTALL_PREFIX}/include/
)
