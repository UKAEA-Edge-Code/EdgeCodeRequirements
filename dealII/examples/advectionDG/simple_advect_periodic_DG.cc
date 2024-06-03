/**
 * DG advection implementation, based on
 * https://github.com/cpraveen/ncmatmw2016/blob/master/2d_scalar_unsteady_legendre/dg.cc
 * .
 */
#include <deal.II/base/function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

constexpr double a_rk[] = {0.0, 3.0 / 4.0, 1.0 / 3.0};
constexpr double b_rk[] = {1.0, 1.0 / 4.0, 2.0 / 3.0};
// Advection velocity
constexpr double vadv_x = 1.0, vadv_y = 0.0;

enum ICsType { gaussian, frankenstein };

/**
 * Convenience function to compute upwind flux
 */
double upwind_flux(double vel, double ul, double ur) {
  if (vel > 0)
    return vel * ul;
  else
    return vel * ur;
}

template <int dim> class InitialCondition : public Function<dim> {
public:
  InitialCondition(ICsType ICs_type) : ICs_type(ICs_type){};
  virtual void
  value_list(const std::vector<Point<dim>> &points, std::vector<double> &values,
             [[maybe_unused]] const unsigned int component = 0) const {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    if (ICs_type == gaussian) {
      constexpr double sigma = 2.0;
      Point<dim> mu = {20.0, 5.0};
      for (unsigned int i = 0; i < values.size(); ++i) {
        values[i] = std::exp(-(points[i] - mu).norm_square() / sigma / sigma);
      }
    } else if (ICs_type == frankenstein) {
      ExcNotImplemented("Frankenstein ICs not implemented yet.");
      for (unsigned int i = 0; i < values.size(); ++i) {
        //     values[i] = 0.0;
      }
    } else {
      AssertThrow(false, ExcMessage("Unknown ICs type"));
    }
  };

private:
  ICsType ICs_type;
};

/**
 * @brief Class to allow boundary condition values to be set explicitly.
 */
template <int dim> class BoundaryValues : public Function<dim> {
public:
  BoundaryValues(){};
  virtual void
  value_list(const std::vector<Point<dim>> &points, std::vector<double> &values,
             [[maybe_unused]] const unsigned int component = 0) const {
    Assert(values.size() == points.size(),
           ExcDimensionMismatch(values.size(), points.size()));

    for (unsigned int i = 0; i < values.size(); ++i) {
      values[i] = 1.0;
    }
  };
};

/**
 * @brief Integrator class using MeshWorker.
 */
template <int dim> class RHSIntegrator {
public:
  RHSIntegrator(const DoFHandler<dim> &dof_handler) : dof_info(dof_handler){};

  MeshWorker::IntegrationInfoBox<dim> int_info;
  MeshWorker::DoFInfo<dim> dof_info;
  MeshWorker::Assembler::ResidualSimple<Vector<double>> assembler;
};

template <int dim> class DGAdvection {
public:
  DGAdvection(unsigned int degree, ICsType ICs_type);
  void run();

private:
  void setup_pBCs();
  void setup_system();
  void assemble_mass_matrix();
  void set_initial_conditions();
  void setup_mesh_worker(RHSIntegrator<dim> &);
  void assemble_rhs(RHSIntegrator<dim> &);
  void compute_dt();
  void solve();
  void output_results(double time);

  Triangulation<dim> triangulation;
  const MappingQ1<dim> mapping;

  unsigned int degree;
  FE_DGP<dim> fe;
  DoFHandler<dim> dof_handler;
  FE_DGP<dim> fe_cell;
  DoFHandler<dim> dof_handler_cell;

  std::vector<Vector<double>> inv_mass_matrix;

  Vector<double> solution;
  Vector<double> solution_old;
  Vector<double> average;
  Vector<double> right_hand_side;
  double dt;
  double cfl;
  ICsType ICs_type;
  double sol_min, sol_max;
  double h_min, h_max;

  std::vector<typename DoFHandler<dim>::cell_iterator> lcell, rcell, bcell,
      tcell;

  typedef MeshWorker::DoFInfo<dim> DoFInfo;
  typedef MeshWorker::IntegrationInfo<dim> CellInfo;

  static void integrate_cell_term(DoFInfo &dinfo, CellInfo &info);
  static void integrate_boundary_term(DoFInfo &dinfo, CellInfo &info);
  static void integrate_face_term(DoFInfo &dinfo1, DoFInfo &dinfo2,
                                  CellInfo &info1, CellInfo &info2);
};

/**
 * Constructor
 */
template <int dim>
DGAdvection<dim>::DGAdvection(unsigned int degree, ICsType ICs_type)
    : mapping(), degree(degree), fe(degree), dof_handler(triangulation),
      fe_cell(0), dof_handler_cell(triangulation), ICs_type(ICs_type) {
  cfl = 0.9 / (2.0 * degree + 1.0);
}

/**
 * @brief Set up periodic boundary conditions in the x direction.
 */
template <int dim> void DGAdvection<dim>::setup_pBCs() {
  // Set up periodicity
  constexpr types::boundary_id lowx_bid = 0, highx_bid = 1;
  constexpr unsigned int xdir_id = 0;
  std::vector<
      GridTools::PeriodicFacePair<typename Triangulation<dim>::cell_iterator>>
      periodic_faces;
  GridTools::collect_periodic_faces(triangulation, lowx_bid, highx_bid, xdir_id,
                                    periodic_faces);
  triangulation.add_periodicity(periodic_faces);
}

template <int dim> void DGAdvection<dim>::setup_system() {

  inv_mass_matrix.resize(triangulation.n_cells());
  for (unsigned int c = 0; c < triangulation.n_cells(); ++c)
    inv_mass_matrix[c].reinit(fe.dofs_per_cell);

  solution.reinit(dof_handler.n_dofs());
  solution_old.reinit(dof_handler.n_dofs());
  right_hand_side.reinit(dof_handler.n_dofs());

  average.reinit(dof_handler_cell.n_dofs());

  assemble_mass_matrix();

  std::cout << "Number of active cells:       "
            << triangulation.n_active_cells() << std::endl;

  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
}

/**
 * Compute and store inverted mass matrix for each cell
 */
template <int dim> void DGAdvection<dim>::assemble_mass_matrix() {
  std::cout << "Constructing mass matrix ..." << std::endl;

  QGauss<dim> quadrature_formula(fe.degree + 1);

  FEValues<dim> fe_values(fe, quadrature_formula,
                          update_values | update_JxW_values);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const unsigned int n_q_points = quadrature_formula.size();

  Vector<double> cell_matrix(dofs_per_cell);

  // Cell iterator
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int c = 0; cell != endc; ++cell, ++c) {
    fe_values.reinit(cell);
    cell_matrix = 0.0;

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        cell_matrix(i) += fe_values.shape_value(i, q_point) *
                          fe_values.shape_value(i, q_point) *
                          fe_values.JxW(q_point);

    // Invert cell_matrix
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      inv_mass_matrix[c](i) = 1.0 / cell_matrix(i);
  }
}

/**
 * Set ICs
 */
template <int dim> void DGAdvection<dim>::set_initial_conditions() {
  InitialCondition<dim> initial_condition(ICs_type);
  VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(fe.degree + 1),
                                      initial_condition, solution);

  // Multiply by inverse mass matrix
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int c = 0; cell != endc; ++cell, ++c) {
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      solution(local_dof_indices[i]) *= inv_mass_matrix[c](i);
  }
}

/**
 * @brief Set up the MeshWorker.
 */
template <int dim>
void DGAdvection<dim>::setup_mesh_worker(RHSIntegrator<dim> &rhs_integrator) {
  std::cout << "Setting up mesh worker ..." << std::endl;

  MeshWorker::IntegrationInfoBox<dim> &int_info = rhs_integrator.int_info;
  MeshWorker::Assembler::ResidualSimple<Vector<double>> &assembler =
      rhs_integrator.assembler;

  const unsigned int n_gauss_points = fe.degree + 1;
  int_info.initialize_gauss_quadrature(n_gauss_points, n_gauss_points,
                                       n_gauss_points);

  // Add solution vector to int_info
  AnyData solution_data;
  solution_data.add<Vector<double> *>(&solution, "solution");
  int_info.cell_selector.add("solution", true, false, false);
  int_info.boundary_selector.add("solution", true, false, false);
  int_info.face_selector.add("solution", true, false, false);

  int_info.initialize_update_flags();
  int_info.add_update_flags_all(update_quadrature_points);
  int_info.add_update_flags_cell(update_gradients);
  int_info.add_update_flags_boundary(update_values);
  int_info.add_update_flags_face(update_values);

  int_info.initialize(fe, mapping, solution_data, Vector<double>());

  // Attach rhs vector to assembler
  AnyData rhs;
  rhs.add<Vector<double> *>(&right_hand_side, "RHS");
  assembler.initialize(rhs);
}

/**
 * Compute min required dt, accounting for spatially varying advection velocity.
 */
template <int dim> void DGAdvection<dim>::compute_dt() {
  std::cout << "Computing time-step ..." << std::endl;

  dt = 1.0e20;
  for (auto cell = dof_handler.begin_active(); cell != dof_handler.end();
       ++cell) {
    double h = cell->diameter() / std::sqrt(2.0);
    double dt_cell = 1.0 / (std::fabs(vadv_x) / h + std::fabs(vadv_y) / h);
    dt = std::min(dt, dt_cell);
  }

  dt *= cfl;
}

template <int dim>
void DGAdvection<dim>::assemble_rhs(RHSIntegrator<dim> &rhs_integrator) {
  right_hand_side = 0.0;

  MeshWorker::loop<dim, dim, MeshWorker::DoFInfo<dim>,
                   MeshWorker::IntegrationInfoBox<dim>>(
      dof_handler.begin_active(), dof_handler.end(), rhs_integrator.dof_info,
      rhs_integrator.int_info, &DGAdvection<dim>::integrate_cell_term,
      &DGAdvection<dim>::integrate_boundary_term,
      &DGAdvection<dim>::integrate_face_term, rhs_integrator.assembler);

  // Multiply by inverse mass matrix
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator cell =
                                                     dof_handler.begin_active(),
                                                 endc = dof_handler.end();
  for (unsigned int c = 0; cell != endc; ++cell, ++c) {
    cell->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      right_hand_side(local_dof_indices[i]) *= inv_mass_matrix[c](i);
  }
}

/**
 * Compute cell integral
 */
template <int dim>
void DGAdvection<dim>::integrate_cell_term(DoFInfo &dinfo, CellInfo &info) {
  const FEValuesBase<dim> &fe_v = info.fe_values();
  const std::vector<double> &sol = info.values[0][0];
  Vector<double> &local_vector = dinfo.vector(0).block(0);
  const std::vector<double> &JxW = fe_v.get_JxW_values();

  Point<dim> vadv = {vadv_x, vadv_y};

  for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) {
    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      local_vector(i) +=
          vadv * fe_v.shape_grad(i, point) * sol[point] * JxW[point];
  }
}

/**
 * Compute boundary integral
 */
template <int dim>
void DGAdvection<dim>::integrate_boundary_term(DoFInfo &dinfo, CellInfo &info) {
  const FEValuesBase<dim> &fe_v = info.fe_values();
  const std::vector<double> &sol = info.values[0][0];

  Vector<double> &local_vector = dinfo.vector(0).block(0);

  const std::vector<double> &JxW = fe_v.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors();

  std::vector<double> g(fe_v.n_quadrature_points);

  static BoundaryValues<dim> boundary_function;
  boundary_function.value_list(fe_v.get_quadrature_points(), g);

  Point<dim> vadv = {vadv_x, vadv_y};
  for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) {
    const double beta_n = vadv * normals[point];
    const double flux = upwind_flux(beta_n, sol[point], g[point]);
    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      local_vector(i) -= flux * fe_v.shape_value(i, point) * JxW[point];
  }
}

/**
 * @brief Integral over internal faces.
 */
template <int dim>
void DGAdvection<dim>::integrate_face_term(DoFInfo &dinfo1, DoFInfo &dinfo2,
                                           CellInfo &info1, CellInfo &info2) {
  const FEValuesBase<dim> &fe_v = info1.fe_values();
  const FEValuesBase<dim> &fe_v_neighbor = info2.fe_values();

  const std::vector<double> &sol1 = info1.values[0][0];
  const std::vector<double> &sol2 = info2.values[0][0];

  Vector<double> &local_vector1 = dinfo1.vector(0).block(0);
  Vector<double> &local_vector2 = dinfo2.vector(0).block(0);

  const std::vector<double> &JxW = fe_v.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors();

  Point<dim> vadv = {vadv_x, vadv_y};
  for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) {
    const double beta_n = vadv * normals[point];
    const double flux = upwind_flux(beta_n, sol1[point], sol2[point]);
    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      local_vector1(i) -= flux * fe_v.shape_value(i, point) * JxW[point];

    for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
      local_vector2(k) +=
          flux * fe_v_neighbor.shape_value(k, point) * JxW[point];
  }
}

/**
 * @brief Set up the MeshWorker, compute dt from the CFL, execute the main time
 * loop, compute error.
 */
template <int dim> void DGAdvection<dim>::solve() {
  RHSIntegrator<dim> rhs_integrator(dof_handler);
  setup_mesh_worker(rhs_integrator);
  compute_dt();

  std::cout << "Start time loop" << std::endl;

  double final_time = 40.0;
  unsigned int iter = 0;
  double time = 0;
  while (time < final_time) {
    // Don't overshoot final_time
    if (time + dt > final_time)
      dt = final_time - time;

    solution_old = solution;

    // 3-stage RK scheme
    for (unsigned int r = 0; r < 3; ++r) {
      assemble_rhs(rhs_integrator);

      for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
        solution(i) = a_rk[r] * solution_old(i) +
                      b_rk[r] * (solution(i) + dt * right_hand_side(i));
    }

    ++iter;
    time += dt;
    std::cout << "Step " << iter << ": t= " << time << std::endl;
    if (std::fmod(iter, 5) == 0 || std::fabs(time - final_time) < 1.0e-14)
      output_results(time);
  }

  // Compare final state to ICs to compute error
  Vector<double> cellwise_errors(triangulation.n_active_cells());
  VectorTools::integrate_difference(
      dof_handler, solution, InitialCondition<dim>(ICs_type), cellwise_errors,
      QGauss<dim>(fe.degree + 2), VectorTools::L2_norm);
  const double u_l2_error = VectorTools::compute_global_error(
      triangulation, cellwise_errors, VectorTools::L2_norm);

  std::cout << "L2 error = " << u_l2_error << std::endl;
}

/**
 * @brief Output current solution to VTK.
 */
template <int dim> void DGAdvection<dim>::output_results(double time) {
  static unsigned int output_num = 0;

  std::string filename =
      "dg-advection-" + Utilities::int_to_string(output_num, 4) + ".vtk";
  std::cout << "Writing solution to <" << filename << ">" << std::endl;
  std::ofstream outfile(filename.c_str());

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "u");
  data_out.build_patches(fe.degree + 1);
  DataOutBase::VtkFlags flags(time, output_num);
  data_out.set_flags(flags);
  data_out.write_vtk(outfile);

  ++output_num;
}

/**
 * Do initialisation and call solve()
 */
template <int dim> void DGAdvection<dim>::run() {
  constexpr unsigned int Nx = 64, Ny = 16;
  constexpr double xmin = 0, xmax = 40, ymin = 0, ymax = 10;
  std::vector<unsigned int> subdivisions = {Nx, Ny};
  const Point<dim> bottom_left = Point<dim>(xmin, ymin);
  const Point<dim> top_right = Point<dim>(xmax, ymax);
  GridGenerator::subdivided_hyper_rectangle(this->triangulation, subdivisions,
                                            bottom_left, top_right, true);

  dof_handler.distribute_dofs(fe);
  dof_handler_cell.distribute_dofs(fe_cell);

  setup_pBCs();
  setup_system();
  set_initial_conditions();

  // Output init state
  output_results(0);

  solve();
}

int main() {
  // Options
  unsigned int degree = 3;
  ICsType ICs_type = gaussian;

  std::cout << "Running with " << MultithreadInfo::n_threads() << " threads"
            << std::endl;

  DGAdvection<2> DG_advection_solver(degree, ICs_type);
  DG_advection_solver.run();

  return 0;
}