add_executable("elasticity-tutorial" elasticity-tutorial.cc)
add_dune_pythonlibs_flags(elasticity-tutorial)
target_link_dune_default_libraries("elasticity-tutorial")

