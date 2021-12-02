#ifndef DUNE_JUPYTER_KERNEL_PRINTGRID_HH
#define DUNE_JUPYTER_KERNEL_PRINTGRID_HH

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/io/file/printgrid.hh>
#include<nlohmann/json.hpp>
#include<xtl/xbase64.hpp>
#include<memory>
#include<functional>

// We inject this into the Dune namespace to enable ADL
namespace Dune {

  template<typename T>
  nlohmann::json mime_bundle_repr(const std::unique_ptr<T>& ptr)
  {
    return mime_bundle_repr(*ptr);
  }

  template<typename T>
  nlohmann::json mime_bundle_repr(const std::shared_ptr<T>& ptr)
  {
    return mime_bundle_repr(*ptr);
  }

  template <class GridView>
  auto mime_bundle_repr(const GridView& gv) -> decltype((gv.grid(), nlohmann::json{})) // currently causes ambiguity
  {
    const auto& grid = gv.grid();
    return mime_bundle_repr(grid);
  }

  // concept check for grid
  template<class Grid>
  auto mime_bundle_repr(const Grid& grid) -> decltype((grid.leafGridView(), nlohmann::json{})) // currently causes ambiguity
  {
    // check for grid dimension
    const int dim = Grid::dimension;

    if constexpr (dim == 2){
      // only print grids with not too many elements
      int gridsize = grid.size(0);
      if (gridsize > 1000){
        auto bundle = nlohmann::json::object();
        bundle["text/plain"] = "Grid is to large for visualization.";
        return bundle;
      }

      // Write out the PNG file using Dune functionality
      int argc;
      char** argv;
      const auto& helper = Dune::MPIHelper::instance(argc, argv);
      Dune::printGrid(grid, helper, "gridexercise");

      // Read it back into a string buffer
      std::ifstream fin("gridexercise_0.png", std::ios::binary);
      std::stringstream buffer;
      buffer << fin.rdbuf();
      fin.close();

      // Return a JSON object for Jupyter
      auto bundle = nlohmann::json::object();
      bundle["image/png"] = xtl::base64encode(buffer.str());
      return bundle;
    }

    auto bundle = nlohmann::json::object();
    bundle["text/plain"] = "Visualization is only available for 2D grids.";
    return bundle;
  }
}

#endif // DUNE_JUPYTER_KERNEL_PRINTGRID_HH

