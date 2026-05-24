#include "args.hpp"
#include <iostream>
#include <sstream>

const char* message=
    "help message, to do";

Args::Args(int argc, char** argv)
  :command(*argv)
{
  char** end = argv + argc;
  ++argv;
  while (argv < end){
    auto option = options.find(*argv);
    if (option != options.end())
      argv = option->second(*this, ++argv, end);
    else{
      std::cerr << "Option " << *argv << " is not valid" << std::endl
        << message << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}

char** Args::print_help(Args& /*_this*/, char** begin, char** /*end*/){
  std::cout << message << std::endl;
  exit(EXIT_SUCCESS);
  return begin;
}

char** Args::set_model(Args& _this, char** begin, char** end){
  if (begin != end)
    _this.model = std::string(*begin);
  else
    abort();
  return ++begin;
}

const std::string& Args::get_model() const{
  return model;
}

char** Args::set_output_file(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.output_file = std::string(*begin);
  else
    abort();
  return ++begin;
}

std::string Args::get_output_file() const{
  return output_file;
}

char** Args::set_delta_t(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.delta_t = std::atof(*begin);
  else
    abort();
  return ++begin;
}

double Args::get_delta_t() const{
  return delta_t;
}

char** Args::set_max_time(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.max_time = std::atof(*begin);
  else
    abort();
  return ++begin;
}

double Args::get_max_time() const{
  return max_time;
}

char** Args::set_mesh_size(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.mesh_size = std::atof(*begin);
  else
    abort();
  return ++begin;
}

double Args::get_mesh_size() const{
  return mesh_size;
}

std::string Args::get_mesh_filename() const{
  if (mesh_size == DEFAULT_MESH_SIZE)
    return DEFAULT_BASE_MESH_FILE + ".msh";
  
  std::string mesh_file_name = DEFAULT_BASE_MESH_FILE + "_h" + (std::ostringstream() << std::fixed << mesh_size).str() + ".msh";

  std::system((std::ostringstream() <<
          "( test ! -f " << mesh_file_name << " || " <<
          "test " << mesh_file_name << " -ot " << DEFAULT_BASE_MESH_FILE + ".geo" << " ) && " << 
          MESH_COMPILING_COMMAND << "-setnumber h " << mesh_size << " -o " << mesh_file_name
        ).str().c_str());

  return mesh_file_name;
}

char** Args::set_theta(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.theta = std::atof(*begin);
  else
    abort();
  return ++begin;
}

double Args::get_theta() const{
  return theta;
}