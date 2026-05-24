#pragma once

#include <unordered_map>
#include <string>
using std::literals::operator""s;

class Args{
  protected:
    const char* command;
    std::unordered_map<std::string, char**(*)(Args&, char**, char**)> options{
      {"-h", print_help},
      {"--help", print_help},
      {"-m", set_model},
      {"--model", set_model},
      {"-o", set_output_file},
      {"--output", set_output_file},
      {"--delta_t", set_delta_t},
      {"--max_time", set_max_time},
      {"--mesh_size", set_mesh_size},
      {"--theta", set_theta}
    };

    std::string model = "TNNP";
    std::string output_file = DEFAULT_BASE_MESH_FILE + ".msh";
    
    double delta_t = 0.005; //ms
    double max_time = 60000; //ms -> 1 simulated hour;
    double mesh_size = DEFAULT_MESH_SIZE; //mm
    double theta = 0.5; //Crank-Nicolson

    static char** print_help(Args& _this, char** begin, char** end);
    static char** set_model(Args& _this, char** begin, char** end);
    static char** set_output_file(Args& _this, char** begin, char** end);
    static char** set_delta_t(Args& _this, char** begin, char** end);
    static char** set_max_time(Args& _this, char** begin, char** end);
    static char** set_mesh_size(Args& _this, char** begin, char** end);
    static char** set_theta(Args& _this, char** begin, char** end);

  public:
    Args(int argc, char** argv);
    const std::string& get_model() const;
    std::string get_output_file() const;
    double get_delta_t() const;
    double get_max_time() const;
    double get_mesh_size() const;
    std::string get_mesh_filename() const;
    double get_theta() const;
};