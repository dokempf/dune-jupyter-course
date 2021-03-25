#ifndef DUNE_JUPYTER_KERNEL_PARAMETERTREE_HH
#define DUNE_JUPYTER_KERNEL_PARAMETERTREE_HH

#include<dune/common/parametertree.hh>
#include<dune/grid/io/file/printgrid.hh>
#include<nlohmann/json.hpp>

// We inject this into the Dune namespace to enable ADL
namespace Dune {

  nlohmann::json mime_bundle_repr(const Dune::ParameterTree& tree)
  {
    std::stringstream stream;
    tree.report(stream);
    auto bundle = nlohmann::json::object();
    bundle["text/plain"] = stream.str();
    return bundle;
  }

}

#endif // DUNE_JUPYTER_KERNEL_PARAMETERTREE_HH
