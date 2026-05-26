#pragma once

#include <unordered_map>
#include <string>
#include <filesystem>
using std::literals::operator""s;

class Args{
  protected:
    const char* command;
    std::unordered_map<std::string, char**(*)(Args&, char**, char**)> options{
      {"-h", print_help},
      {"--help", print_help},
      {"-m", set_model},
      {"--model", set_model},
      {"-o", set_output_file_name},
      {"--output", set_output_file_name},
      {"--delta_t", set_delta_t},
      {"--max_time", set_max_time},
      {"--mesh_size", set_mesh_size},
      {"--theta", set_theta},
      {"--implicit_euler", set_implicit_euler},
      {"--explicit_euler", set_explicit_euler},
      {"--crank_nicolson", set_crank_nicolson}
    };

    std::string model = DEFAULT_MODEL;
    std::string output_file_name = DEFAULT_OUTPUT_FILE;
    
    double delta_t = DEFAULT_DELTA_T; //ms
    double max_time = DEFAULT_MAX_TIME; //ms;
    double mesh_size = DEFAULT_MESH_SIZE; //mm
    double theta = DEFAULT_THETA; 

    void print_help_and_exit(int exit_status);

    static char** print_help(Args& _this, char** begin, char** end);
    static char** set_model(Args& _this, char** begin, char** end);
    static char** set_output_file_name(Args& _this, char** begin, char** end);
    static char** set_delta_t(Args& _this, char** begin, char** end);
    static char** set_max_time(Args& _this, char** begin, char** end);
    static char** set_mesh_size(Args& _this, char** begin, char** end);
    static char** set_theta(Args& _this, char** begin, char** end);
    static char** set_implicit_euler(Args& _this, char** begin, char** end);
    static char** set_explicit_euler(Args& _this, char** begin, char** end);
    static char** set_crank_nicolson(Args& _this, char** begin, char** end);


  public:
    Args(int argc, char** argv);
    const std::string& get_model() const;
    std::string get_output_file_name() const;
    double get_delta_t() const;
    double get_max_time() const;
    double get_mesh_size() const;
    std::string get_mesh_filename() const;
    double get_theta() const;

    friend std::ostream& operator<<(std::ostream& os, const Args& args);
};