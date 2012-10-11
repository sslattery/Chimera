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
  Chimera_Assertion.hpp
  Chimera_BoostRNG.hpp
  Chimera_History.hpp
  Chimera_History_def.hpp
  Chimera_JacobiSplit.hpp
  Chimera_JacobiSplit_def.hpp
  Chimera_LinearOperatorSplit.hpp
  Chimera_LinearOperatorSplitFactory.hpp
  Chimera_LinearOperatorSplitFactory_def.hpp
  Chimera_LinearProblem.hpp
  Chimera_LinearProblem_def.hpp
  Chimera_RNGTraits.hpp
  Chimera_StationaryIteration.hpp
  Chimera_StationaryIteration_def.hpp
  ) 

APPEND_SET(SOURCES
  Chimera_Assertion.cpp
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
