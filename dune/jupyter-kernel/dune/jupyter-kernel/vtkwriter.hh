#ifndef DUNE_JUPYTER_KERNEL_VTKWRITER_HH
#define DUNE_JUPYTER_KERNEL_VTKWRITER_HH

#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<nlohmann/json.hpp>


namespace Dune {

  template<typename GV>
  nlohmann::json mime_bundle_repr(Dune::VTKWriter<GV>& vtkwriter)
  {
    // Write out the VTK file into a temporary file
    vtkwriter.write("output", Dune::VTK::appendedraw);
    auto bundle = nlohmann::json::object();
    bundle["text/html"] = "VTK Output has been written to file. You can <a href=\"output.vtu\">download it</a> and visualize it in ParaView.";
    return bundle;
  }

  template<typename GV>
  nlohmann::json mime_bundle_repr(Dune::SubsamplingVTKWriter<GV>& vtkwriter)
  {
    // Write out the VTK file into a temporary file
    vtkwriter.write("output", Dune::VTK::appendedraw);
    auto bundle = nlohmann::json::object();
    bundle["text/html"] = "VTK Output has been written to file. You can <a href=\"output.vtu\">download it</a> and visualize it in ParaView.";
    return bundle;
  }


} // namespace Dune

#endif