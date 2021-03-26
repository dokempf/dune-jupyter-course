Kernel stuff:
* Update to new version of xeus-cling when out
* Try if `std::call_once` issue is fixed on new version

CMake integration stuff:
* Carry over compile flags from the Dune buildsystem (esp. defines)
* Check why we cannot move `dune/jupyter.hh` into the `xeus_cling.hh` header

Dune/Jupyter integration stuff
* Add a concept check around the grid print so that it only prints grids
* Check `dim == 2` before doing `Dune::printgrid`
* Limit the Gnuplot printing thing to a certain cells threshold, beyond that VTK needs to be used
* Check links for Paraview integration:
  * https://www.paraview.org/pipermail/paraview/2017-January/038984.html

Jupyter Stuff:
* Check exercise extension
* Check why precommit hook sometime fails
