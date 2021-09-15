#ifndef DUNE_JUPYTER_KERNEL_VTKWRITER_HH
#define DUNE_JUPYTER_KERNEL_VTKWRITER_HH

#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<nlohmann/json.hpp>


namespace Dune {

  template<typename GV>
  nlohmann::json mime_bundle_repr(Dune::VTKWriter<GV>& vtkwriter)
  {
    //read name from temporary file             // this can be created using a magic command "%%file filename" in a cell, rest of cell content is then written in this file
    std::string outputname;
    std::ifstream myfile("name.txt");           // so in order for this to work the magic command "%%file name.txt" is necessary, followed by the prefered output name
    if (myfile.is_open())
    {
      std::getline(myfile,outputname);
      myfile.close();
    }
    else
      outputname = "output";                    // if no magic command is provided, it still works as before

    // Write out the VTK file into a temporary file
    vtkwriter.write(outputname, Dune::VTK::appendedraw);
    auto bundle = nlohmann::json::object();
    std::string bundletxt = "VTK Output has been written to file. You can <a href=\""+ outputname + ".vtu\">download it</a> and visualize it in ParaView.";
    bundle["text/html"] = bundletxt;
    return bundle;
  }

  template<typename GV>
  nlohmann::json mime_bundle_repr(Dune::SubsamplingVTKWriter<GV>& vtkwriter)
  {
    //read name from temporary file
    std::string outputname;
    std::ifstream myfile("name.txt");
    if (myfile.is_open())
    {
      std::getline(myfile,outputname);
      myfile.close();
    }
    else
      outputname = "output";

    // Write out the VTK file into a temporary file
    vtkwriter.write(outputname, Dune::VTK::appendedraw);
    auto bundle = nlohmann::json::object();
    std::string bundletxt = "VTK Output has been written to file. You can <a href=\""+ outputname + ".vtu\">download it</a> and visualize it in ParaView.";
    bundle["text/html"] = bundletxt;
    return bundle;
  }

  template<typename GV>
  nlohmann::json mime_bundle_repr(Dune::VTKSequenceWriter<GV>& vtkwriter)
  {
    //read name from temporary file
    std::string outputname;
    std::ifstream myfile("name.txt");
    if (myfile.is_open())
    {
      std::getline(myfile,outputname);
      myfile.close();
    }
    else
      outputname = "output";

    //use system commands to zip folder
    std::string cmd = "zip -r " + outputname + ".zip " + outputname;
    std::system(cmd.c_str());

    auto bundle = nlohmann::json::object();
    std::string bundletxt = "VTK Output has been written to folder. You can <a href=\""+ outputname+ ".zip\"> download it </a> and the <a href=\""+ outputname + ".pvd\"> corresponding file</a> .";
    bundle["text/html"] = bundletxt;
    return bundle;
  }

} // namespace Dune

#endif