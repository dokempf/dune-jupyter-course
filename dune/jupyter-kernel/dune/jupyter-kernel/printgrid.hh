#ifndef DUNE_JUPYTER_KERNEL_PRINTGRID_HH
#define DUNE_JUPYTER_KERNEL_PRINTGRID_HH

#include<dune/common/parallel/mpihelper.hh>
#include<dune/grid/io/file/printgrid.hh>
#include<nlohmann/json.hpp>
#include<xtl/xbase64.hpp>

// We inject this into the Dune namespace to enable ADL
namespace Dune {

  // TODO: Add a concept check for a Dune grid
  template<typename Grid>
  nlohmann::json mime_bundle_repr(const std::unique_ptr<Grid>& grid)
  {
    //check for grid dimension
    const int dim = Grid::dimension;
    if(dim > 2){
      auto bundle = nlohmann::json::object();
      bundle["text/plain"] = "Visualization is only available for 2D grids.";
      return bundle;
    }
    // Write out the PNG file using Dune functionality
    int argc;
    char** argv;
    const auto& helper = Dune::MPIHelper::instance(argc, argv);
    Dune::printGrid(*grid, helper, "gridexercise");

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

/*
  template<int dim, typename T, typename U, template<int,typename> typename Grid>
  nlohmann::json mime_bundle_repr(const std::unique_ptr<Grid<dim,T>, U > & grid)
  {
    // Write out the PNG file using Dune functionality
    int argc;
    char** argv;
    const auto& helper = Dune::MPIHelper::instance(argc, argv);
    Dune::printGrid(*grid, helper, "gridexercise");

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

  //if grid has dim=3 it cannot be printed
  template<typename T, typename U, template<int,typename> typename Grid>
  nlohmann::json mime_bundle_repr(const std::unique_ptr<Grid<3,T>, U > & grid)
  {
      auto bundle = nlohmann::json::object();
      bundle["text/html"] = "Visualization is only available for 2D grids.";
      return bundle;
  }
*/

  //the function grid() of a gridview returns a reference to a grid, but this should be done using concepts
  /*template<int dim, typename T, template<int,typename> typename Grid>
  nlohmann::json mime_bundle_repr(const Grid<dim,T>& grid)
  {
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

  template<typename T, template<int,typename> typename Grid>
  nlohmann::json mime_bundle_repr(const Grid<3,T>& grid)
  {
      auto bundle = nlohmann::json::object();
      bundle["text/html"] = "Visualization is only available for 2D grids.";
      return bundle;
  }

  //printgrid function for grid
  template <typename Grid>
  nlohmann::json mime_bundle_repr(const GridView<Grid> & gridview)
  {
      const auto& grid = gridview.grid();
      return mime_bundle_repr(grid);
  }
  */
}

#endif // DUNE_JUPYTER_KERNEL_PRINTGRID_HH

