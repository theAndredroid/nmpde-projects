#include "Current.hpp"
struct MVParameters
{
  // --- Valori limite del voltaggio ---
  double u_o        = 0.0;    // voltaggio a riposo (adimensionale)
  double u_u        = 1.61;   // voltaggio massimo upstroke

  // --- Soglie per le funzioni di Heaviside ---
  double theta_v    = 0.3;    // soglia per J_fi e gate v
  double theta_w    = 0.13;   // soglia per J_so, J_si, gate w
  double theta_v_m  = 0.1;    // soglia per tau_v- (quale ramo)
  double theta_o    = 0.005;  // soglia per tau_o e w_inf
 
  // --- Parametri per tau_v- (costante di tempo v in chiusura) ---
  double tau_v1m    = 80.0;   // tau_v - quando u < theta_vm
  double tau_v2m    = 1.4506; // tau_v- quando u > theta_vm
 
  // --- Parametri per tau_v+ (costante di tempo v in apertura) ---
  double tau_v_p    = 1.4506;
 
  // --- Parametri per tau_w- (costante di tempo w in chiusura) ---
  double tau_w1_m   = 70.0;   // valore minimo di tau_w-
  double tau_w2_m   = 8.0;    // valore massimo di tau_w-
  double k_w_m      = 200.0;  // slope della sigmoide per tau_w-
  double u_w_m      = 0.016;  // punto di mezzo della sigmoide
 
  // --- Parametri per tau_w+ (costante di tempo w in apertura) ---
  double tau_w_p    = 280.0;
 
  // --- Parametri per J_fi (fast inward - sodio) ---
  double tau_fi     = 0.078;
 
  // --- Parametri per J_so (slow outward - potassio) ---
  double tau_o1     = 410.0;  // tau_o quando u < theta_o
  double tau_o2     = 7.0;    // tau_o quando u > theta_o
  double tau_so1    = 91.0;   // tau_so minimo
  double tau_so2    = 0.8;    // tau_so massimo
  double k_so       = 2.1;    // slope della sigmoide per tau_so
  double u_so       = 0.6;    // punto di mezzo della sigmoide
 
  // --- Parametri per s (quarta variabile, morfologia AP) ---
  double tau_s1     = 2.7342; // tau_s quando u < theta_w
  double tau_s2     = 4.0;    // tau_s quando u > theta_w
  double k_s        = 2.0994; // slope della tanh per s_inf
  double u_s        = 0.9087; // punto di mezzo della tanh
 
  // --- Parametri per J_si (slow inward - calcio) ---
  double tau_si     = 3.3849;
 
  // --- Parametri per w_inf ---
  double tau_w_inf  = 0.01;   // usato nel calcolo di w_inf
  double w_inf_star = 0.5;    // valore di w_inf quando u > theta_o
 
  // --- Diffusione ---
  // D = 1.171 cm^2/s dal paper (Appendice A)
  // Qui in unità adimensionali del modello
  double D = 1.171e-4;        // adattato alle unità del problema
};
 
// ============================================================

void
Current::setup()
{
  pcout << "===============================================" << std::endl;

  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    // Read serial mesh.
    Triangulation<dim> mesh_serial;

    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(mesh_serial);

      std::ifstream mesh_file(mesh_file_name);
      grid_in.read_msh(mesh_file);
    }

    // Copy the serial mesh into the parallel one.
    {
      GridTools::partition_triangulation(mpi_size, mesh_serial);

      const auto construction_data = TriangulationDescription::Utilities::
        create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
      mesh.create_triangulation(construction_data);
    }

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_Q<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGauss<dim>>(r + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
    const IndexSet locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);

    pcout << "  Initializing the sparsity pattern" << std::endl;
    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    pcout << "  Initializing the system matrix" << std::endl;
    system_matrix.reinit(sparsity);

    pcout << "  Initializing vectors" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    v_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    w_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    s_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    // v.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    // w.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    // s.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  }
}

// void Current::compute_ionic_currents(){
  
// }

// void
// Current::solve_ionic()
// {
//   MVParameters p; // parametri del modello MV

//   // Iteriamo su tutti i DoF locali
//   const IndexSet &locally_owned = dof_handler.locally_owned_dofs();

//   for (const auto i : locally_owned)
//     {
//       // Leggi u, v, w, s al timestep precedente
//       const double u = solution_old(i);
//       const double v = var_v(i);
//       const double w = var_w(i);
//       const double s = var_s(i);

//       // -------------------------------------------------------
//       // Funzioni di Heaviside H(u - soglia)
//       // H = 1 se u > soglia, 0 altrimenti
//       // Governano quali canali sono aperti o chiusi
//       // -------------------------------------------------------
//       const double Hv  = (u > p.theta_v)  ? 1.0 : 0.0;
//       const double Hw  = (u > p.theta_w)  ? 1.0 : 0.0;
//       const double Hvm = (u > p.theta_v_m) ? 1.0 : 0.0;
//       const double Ho  = (u > p.theta_o)  ? 1.0 : 0.0;

//       // -------------------------------------------------------
//       // Costanti di tempo dipendenti da u
//       // Cambiano in base allo stato del voltaggio
//       // -------------------------------------------------------

//       // tau_v-: controlla velocita' di chiusura del canale Na
//       const double tau_vm = (1.0 - Hvm) * p.tau_v1m
//                           +        Hvm  * p.tau_v2m;

//       // tau_w-: sigmoide, controlla chiusura del canale Ca
//       const double tau_wm = p.tau_w1_m
//                           + (p.tau_w2_m - p.tau_w1_m)
//                           * (1.0 + std::tanh(p.k_w_m * (u - p.u_w))) / 2.0;

//       // tau_s: controlla la quarta variabile s
//       const double tau_s = (1.0 - Hw) * p.tau_s1
//                          +        Hw  * p.tau_s2;

//       // -------------------------------------------------------
//       // Valori asintotici (a cosa tendono v, w, s)
//       // -------------------------------------------------------

//       // v_inf: a riposo (u basso) v tende a 1 (canale pronto)
//       //        eccitato (u alto) v tende a 0 (canale inattivato)
//       const double v_inf = (u > p.theta_vm) ? 0.0 : 1.0;

//       // w_inf: dipende da u tramite theta_o
//       const double w_inf = (1.0 - Ho) * (1.0 - u / p.tau_winf)
//                          +        Ho  * p.w_inf_star;

//       // s_inf: sigmoide centrata in u_s
//       const double s_inf = (1.0 + std::tanh(p.k_s * (u - p.u_s))) / 2.0;

//       // -------------------------------------------------------
//       // Euler esplicito per le tre ODE
//       //
//       // dv/dt = (1-Hv)*(v_inf - v)/tau_vm - Hv*v/tau_vp
//       //   - quando u < theta_v: v recupera verso v_inf
//       //   - quando u > theta_v: v decade (canale si inattiva)
//       //
//       // dw/dt = (1-Hw)*(w_inf - w)/tau_wm - Hw*w/tau_wp
//       //   - stessa logica per il canale calcio
//       //
//       // ds/dt = (s_inf - s) / tau_s
//       //   - s segue sempre s_inf con costante tau_s
//       // -------------------------------------------------------
//       var_v(i) += delta_t * ((1.0 - Hv) * (v_inf - v) / tau_vm
//                             -        Hv  *  v          / p.tau_vp);

//       var_w(i) += delta_t * ((1.0 - Hw) * (w_inf - w) / tau_wm
//                             -        Hw  *  w          / p.tau_wp);

//       var_s(i) += delta_t * (s_inf - s) / tau_s;
//     }

//   // Comunica i valori aggiornati tra i processi MPI
//   var_v.compress(VectorOperation::insert);
//   var_w.compress(VectorOperation::insert);
//   var_s.compress(VectorOperation::insert);
// }

// ============================================================

// void
// Current::assemble()
// {
//   // Number of local DoFs for each element.
//   const unsigned int dofs_per_cell = fe->dofs_per_cell;

//   // Number of quadrature points for each element.
//   const unsigned int n_q = quadrature->size();

//   FEValues<dim> fe_values(*fe,
//                           *quadrature,
//                           update_values | update_gradients |
//                             update_quadrature_points | update_JxW_values);

//   // Local matrix and vector.
//   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
//   Vector<double>     cell_rhs(dofs_per_cell);

//   std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

//   // Reset the global matrix and vector, just in case.
//   system_matrix = 0.0;
//   system_rhs    = 0.0;

//   // Evaluation of the old solution on quadrature nodes of current cell.
//   std::vector<double> solution_old_values(n_q);

//   // Evaluation of the gradient of the old solution on quadrature nodes of
//   // current cell.
//   std::vector<Tensor<1, dim>> solution_old_grads(n_q);

//   for (const auto &cell : dof_handler.active_cell_iterators())
//     {
//       if (!cell->is_locally_owned())
//         continue;

//       fe_values.reinit(cell);

//       cell_matrix = 0.0;
//       cell_rhs    = 0.0;

//       // Evaluate the old solution and its gradient on quadrature nodes.
//       fe_values.get_function_values(solution, solution_old_values);
//       fe_values.get_function_gradients(solution, solution_old_grads);

//       for (unsigned int q = 0; q < n_q; ++q)
//         {
//           const double mu_loc = mu(fe_values.quadrature_point(q));

//           const double f_old_loc =
//             f(fe_values.quadrature_point(q), time - delta_t);
//           const double f_new_loc = f(fe_values.quadrature_point(q), time);

//           for (unsigned int i = 0; i < dofs_per_cell; ++i)
//             {
//               for (unsigned int j = 0; j < dofs_per_cell; ++j)
//                 {
//                   // Time derivative.
//                   cell_matrix(i, j) += (1.0 / delta_t) *             //
//                                        fe_values.shape_value(i, q) * //
//                                        fe_values.shape_value(j, q) * //
//                                        fe_values.JxW(q);

//                   // Diffusion.
//                   cell_matrix(i, j) +=
//                     theta * mu_loc *                             //
//                     scalar_product(fe_values.shape_grad(i, q),   //
//                                    fe_values.shape_grad(j, q)) * //
//                     fe_values.JxW(q);
//                 }

//               // Time derivative.
//               cell_rhs(i) += (1.0 / delta_t) *             //
//                              fe_values.shape_value(i, q) * //
//                              solution_old_values[q] *      //
//                              fe_values.JxW(q);

//               // Diffusion.
//               cell_rhs(i) -= (1.0 - theta) * mu_loc *                   //
//                              scalar_product(fe_values.shape_grad(i, q), //
//                                             solution_old_grads[q]) *    //
//                              fe_values.JxW(q);

//               // Forcing term.
//               cell_rhs(i) +=
//                 (theta * f_new_loc + (1.0 - theta) * f_old_loc) * //
//                 fe_values.shape_value(i, q) *                     //
//                 fe_values.JxW(q);
//             }
//         }

//       cell->get_dof_indices(dof_indices);

//       system_matrix.add(dof_indices, cell_matrix);
//       system_rhs.add(dof_indices, cell_rhs);
//     }

//   system_matrix.compress(VectorOperation::add);
//   system_rhs.compress(VectorOperation::add);

//   // Homogeneous Neumann boundary conditions: we do nothing.
// }

void
Current::solve_linear_system()
{
  TrilinosWrappers::PreconditionSSOR preconditioner;
  preconditioner.initialize(
    system_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  ReductionControl solver_control(/* maxiter = */ 10000,
                                  /* tolerance = */ 1.0e-16,
                                  /* reduce = */ 1.0e-6);

  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);

  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  pcout << solver_control.last_step() << "GMRES iterations" << std::endl;
}

void
Current::output() const
{
  DataOut<dim> data_out;

  data_out.add_data_vector(dof_handler, solution, "membrane tension");
  data_out.add_data_vector(dof_handler, v_owned, "v");
  data_out.add_data_vector(dof_handler, w_owned, "w");
  data_out.add_data_vector(dof_handler, s_owned, "s");

  // Add vector for parallel partition.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string output_file_name = "output-" + mesh_path.stem().string();

  data_out.write_vtu_with_pvtu_record(/* folder = */ "./",
                                      /* basename = */ output_file_name,
                                      /* index = */ timestep_number,
                                      MPI_COMM_WORLD);
}

void
Current::run()
{
  // Setup initial conditions.
  {
    setup();

    VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), solution_owned);
    solution = solution_owned;

    // -------------------------------------------------------
    // Iniziliaze the auxilial values to initial conditions.
    //   v = 1 -> canale sodium channel ready to be opened
    //   w = 1 -> calcium channel ready to be opened  
    //   s = 0 -> 4-th setted to zero
    // -------------------------------------------------------
    VectorTools::interpolate(dof_handler, Functions::ConstantFunction<dim>(1.0), v_owned);
    VectorTools::interpolate(dof_handler, Functions::ConstantFunction<dim>(1.0), w_owned);
    VectorTools::interpolate(dof_handler, Functions::ZeroFunction<dim>(), s_owned);

    pcout << "  Ionic variables initialized (v=1, w=1, s=0)" << std::endl;

    // compute_ionic_currents();

    // time            = 0.0;
    // timestep_number = 0;

    // Output initial condition.
    output();
  }

  pcout << "===============================================" << std::endl;

  // Time-stepping loop.
  // while (time < T - 0.5 * delta_t)
  //   {
  //     time += delta_t;
  //     ++timestep_number;

  //     pcout << "Timestep " << std::setw(3) << timestep_number
  //           << ", time = " << std::setw(4) << std::fixed << std::setprecision(2)
  //           << time << " : ";

  //     assemble();
  //     solve_linear_system();

  //     // Perform parallel communication to update the ghost values of the
  //     // solution vector.
  //     solution = solution_owned;

  //     output();
  //   }
}