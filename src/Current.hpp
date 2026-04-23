#ifndef CURRENT_HPP
#define CURRENT_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// #include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace dealii;

/**
 * Class managing the differential problem.
 */
class Current
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Constructor.
  Current(const std::string                            &mesh_file_name_,
       const unsigned int                              &r_,
       const double                                    &T_,
       const double                                    &theta_,
       const double                                    &delta_t_)
    : mesh_file_name(mesh_file_name_)
    , r(r_)
    , T(T_)
    , theta(theta_)
    , delta_t(delta_t_)
    , mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , mesh(MPI_COMM_WORLD)
    , pcout(std::cout, mpi_rank == 0)
  {}

  // Run the time-dependent simulation.
  void
  run();

protected:
  // Initialization.
  void
  setup();

  // System assembly.
  void
  assemble();

  // System solution.
  void
  solve_linear_system();

  void integrate_auxiliar_variables();

  void compute_ionic_currents();

  void activation_time();

  // Output.
  void
  output() const;

  void
  output_activation_time() const;

  // Name of the mesh.
  const std::string mesh_file_name;

  // Polynomial degree.
  const unsigned int r;

  // Final time.
  const double T;

  // Theta parameter for the theta method.
  const double theta;

  // Time step.
  const double delta_t;

  // Current time in ms.
  double time = 0.0;

  // Current timestep number.
  unsigned int timestep_number = 0;

  // // Forcing term.
  // std::function<double(const Point<dim> &, const double &)> f;

  // Number of MPI processes.
  const unsigned int mpi_size;

  // Rank of the current MPI process.
  const unsigned int mpi_rank;

  // Triangulation.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // System matrix.
  TrilinosWrappers::SparseMatrix system_matrix;

  // System right-hand side.
  TrilinosWrappers::MPI::Vector system_rhs;

  // System solution, without ghost elements.
  TrilinosWrappers::MPI::Vector solution_owned;

  // System solution, with ghost elements.
  TrilinosWrappers::MPI::Vector solution;

  //
  TrilinosWrappers::MPI::Vector Time;

  // Currents, without ghost elements.
  TrilinosWrappers::MPI::Vector J_fi_owned;
  TrilinosWrappers::MPI::Vector J_so_owned;
  TrilinosWrappers::MPI::Vector J_si_owned;

  // Currents, with ghost elements.
  TrilinosWrappers::MPI::Vector J_fi;
  TrilinosWrappers::MPI::Vector J_so;
  TrilinosWrappers::MPI::Vector J_si;

  // Auxiliar variables, without ghost elements.
  TrilinosWrappers::MPI::Vector v_owned;
  TrilinosWrappers::MPI::Vector w_owned;
  TrilinosWrappers::MPI::Vector s_owned;

  // Auxiliar variables, with ghost elements.
  // TrilinosWrappers::MPI::Vector v;
  // TrilinosWrappers::MPI::Vector w;
  // TrilinosWrappers::MPI::Vector s;

  // Output stream for process 0.
  ConditionalOStream pcout;
};

#endif