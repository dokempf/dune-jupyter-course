#ifndef DUNE_JUPYTER_KERNEL_PRINTREFELEM_HH
#define DUNE_JUPYTER_KERNEL_PRINTREFELEM_HH

#include<nlohmann/json.hpp>

// We inject this into the Dune::Geo namespace to enable ADL
namespace Dune::Geo {
  template <typename T, int dim>
  nlohmann::json mime_bundle_repr(const Dune::Geo::ReferenceElement<Dune::Geo::ReferenceElementImplementation<T, dim>>& refelem)
  {
    //get the type of the refelem
    const Dune::GeometryType& geotype = refelem.type();

    //get the matching filename for the geometry type
    std::string filename = "referenceElements/";
    if(geotype.isLine())
      filename += "gg_line.svg";
    else if(geotype.isTriangle())
      filename += "gg_triangle.svg";
    else if(geotype.isQuadrilateral())
      filename += "gg_quadrilateral.svg";
    else if(geotype.isTetrahedron())
      filename += "gg_tetrahedron.svg";
    else if(geotype.isPyramid())
      filename += "gg_pyramid.svg";
    else if(geotype.isPrism())
      filename += "gg_prism.svg";
    else if(geotype.isHexahedron())
      filename += "gg_hexahedron.svg";
    else {
      //if it is non of the above, return message
      auto bundle = nlohmann::json::object();
      bundle["text/html"] = "No visualization available for this reference element.<br> Available options are: line, triangle, quadrilateral, tetrahedron, pyramid, prism and hexahedron.";
      return bundle;
    }

    std::ifstream fin(filename, std::ios::binary);
    std::stringstream buffer;
    buffer << fin.rdbuf();
    fin.close();

    // Return a JSON object for Jupyter
    auto bundle = nlohmann::json::object();
    bundle["image/svg+xml"] = buffer.str(); //.svg does not need to be encoded in base64
    return bundle;
  }
}

#endif // DUNE_JUPYTER_KERNEL_PRINTREFELEM_HH