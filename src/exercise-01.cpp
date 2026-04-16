#include "Current.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  // constexpr unsigned int dim = Current::dim;

  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);

  Current problem(/*mesh_filename = */ "mesh/mesh-square.msh",
               /* degree = */ 1,
               /* T = */ 3.0,
               /* theta = */ 0.5, //crank-Nicolson
               /* delta_t = */ 0.05);

  problem.run();

  return 0;
}