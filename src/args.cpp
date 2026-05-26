#include "args.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

#include <deal.II/base/mpi.h>

Args::Args(int argc, char** argv)
  :command(*argv)
{
  char** end = argv + argc;
  ++argv;
  while (argv < end){
    auto option = options.find(*argv);
    if (option != options.end())
      argv = option->second(*this, ++argv, end);
    else
      this->print_help_and_exit(EXIT_FAILURE);
  }
}

void Args::print_help_and_exit(int exit_status)
{
  const unsigned int mpi_rank = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (mpi_rank == 0){
    std::ostream &os = (exit_status == EXIT_SUCCESS)? std::cout : std::cerr;
    
    std::string available_models = "";
    std::ifstream model_file("models/models.csv");
    if (model_file.is_open()) {
      std::string line;
      std::getline(model_file, line); // skip header
      while (std::getline(model_file, line)) {
        std::stringstream ss(line);
        std::string model_name;
        if (std::getline(ss, model_name, ',')) {
          // Trim whitespace
          model_name.erase(0, model_name.find_first_not_of(" \t"));
          model_name.erase(model_name.find_last_not_of(" \t") + 1);
          if (!available_models.empty()) available_models += ", ";
          available_models += model_name;
        }
      }
    }

    os << "Usage: " << command << " [options]\n"
      << "Options:\n"
      << "  -h, --help            Show this help message and exit\n"
      << "  -m, --model <name>    Set the model name (default: " << model << ")\n";
    if (!available_models.empty()) {
      os << "                        Available models: " << available_models << "\n";
    }
    os << "  -o, --output <file>   Set the output file name\n"
      << "  --delta_t <val>       Set delta_t in ms (default: " << delta_t << ")\n"
      << "  --max_time <val>      Set max_time in ms (default: " << max_time << ")\n"
      << "  --mesh_size <val>     Set mesh_size in mm (default: " << mesh_size << ")\n"
      << "  --theta <val>         Set theta for theta-method for time discretization (default: " << theta << ")\n"
      << "  --implicit_euler      Set theta = 1.0\n"
      << "  --explicit_euler      Set theta = 0.0\n"
      << "  --crank_nicolson      Set theta = 0.5\n";
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  std::exit(exit_status);
}

char **Args::print_help(Args & _this, char **begin, char ** /*end*/)
{
  _this.print_help_and_exit(EXIT_SUCCESS);
  return begin;
}

char** Args::set_model(Args& _this, char** begin, char** end){
  if (begin != end)
    _this.model = std::string(*begin);
  else
    _this.print_help_and_exit(EXIT_FAILURE);
  return ++begin;
}

const std::string& Args::get_model() const{
  return model;
}

char** Args::set_output_file_name(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.output_file_name = std::string(*begin);
  else
    _this.print_help_and_exit(EXIT_FAILURE);
  return ++begin;
}

std::string Args::get_output_file_name() const{
  return output_file_name;
}

char** Args::set_delta_t(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.delta_t = std::atof(*begin);
  else
    _this.print_help_and_exit(EXIT_FAILURE);
  return ++begin;
}

double Args::get_delta_t() const{
  return delta_t;
}

char** Args::set_max_time(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.max_time = std::atof(*begin);
  else
    _this.print_help_and_exit(EXIT_FAILURE);
  return ++begin;
}

double Args::get_max_time() const{
  return max_time;
}

char** Args::set_mesh_size(Args &_this, char **begin, char **end){
  if (begin != end)
    _this.mesh_size = std::atof(*begin);
  else
    _this.print_help_and_exit(EXIT_FAILURE);
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
    _this.print_help_and_exit(EXIT_FAILURE);
  return ++begin;
}

char **Args::set_implicit_euler(Args &_this, char **begin, char **end)
{
  _this.theta = 1.0;
  return begin;
}

char **Args::set_explicit_euler(Args &_this, char **begin, char **end)
{
  _this.theta = 0.0;
  return begin;
}

char **Args::set_crank_nicolson(Args &_this, char **begin, char **end)
{
  _this.theta = 0.5;
  return begin;
}

double Args::get_theta() const{
  return theta;
}

std::ostream& operator<<(std::ostream& os, const Args& args){
  os << "Configuration:" << std::endl
     << "  Model:       " << args.model << std::endl
     << "  Output file: " << args.output_file_name << std::endl
     << "  Delta t:     " << args.delta_t << " ms" << std::endl
     << "  Max time:    " << args.max_time << " ms" << std::endl
     << "  Mesh size:   " << args.mesh_size << " mm" << std::endl
     << "  Theta:       " << args.theta << std::endl;
  return os;
}