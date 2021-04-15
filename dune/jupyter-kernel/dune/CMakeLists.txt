add_subdirectory(jupyter-kernel)

install(
  FILES
    common.hh
    geometry.hh
    grid.hh
    jupyter.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune
)
