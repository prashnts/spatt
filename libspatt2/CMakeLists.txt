SET(spatt_lib_src alphabet.cc dfa.cc gstat.cc markov.cc nfa.cc pmc.cc sequence.cc stat.cc vmarkov.cc xstat.cc xwaiting.cc)

ADD_LIBRARY(spatt-lib-static STATIC ${spatt_lib_src})
SET_TARGET_PROPERTIES(spatt-lib-static PROPERTIES OUTPUT_NAME spatt2)

OPTION(WITH_SHARED_LIB "Compile shared library" OFF)
IF(WITH_SHARED_LIB)
  MESSAGE(STATUS "Compilation of the shared library activated")
  ADD_LIBRARY(spatt-lib SHARED ${spatt_lib_src} spattc.cc)
  TARGET_LINK_LIBRARIES(spatt-lib m)
  SET_TARGET_PROPERTIES(spatt-lib PROPERTIES OUTPUT_NAME spatt2)
  INSTALL(TARGETS spatt-lib DESTINATION "lib")
ENDIF(WITH_SHARED_LIB)
