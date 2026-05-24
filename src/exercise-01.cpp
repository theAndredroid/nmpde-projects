#include "Current.hpp"
#include "args.hpp"

// Main function.
int main(int argc, char *argv[]){
  // constexpr unsigned int dim = Current::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  Args args(argc, argv);
  
  std::cout << args.get_mesh_filename() << std::endl;

  Current problem(/*mesh_filename = */ args.get_mesh_filename(),
               /* degree = */ 1,
               /* T = */ args.get_max_time(),
               /* theta = */ args.get_theta(),
               /* delta_t = */ args.get_delta_t());
  
  problem.set_model(args.get_model());

  problem.run();

  return 0;
}
