#install headers
add_subdirectory("materials")


install(FILES elasticity-tutorial.hh 
        linearElasticityDisplacementAssembler.hh
        meshUtilities.hh
        statistics.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/elasticity-tutorial)