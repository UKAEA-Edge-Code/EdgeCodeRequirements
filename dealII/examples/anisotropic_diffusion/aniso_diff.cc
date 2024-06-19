#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <filesystem>
#include <math.h>
namespace fs = std::filesystem;
using namespace dealii;

const fs::path this_dir = fs::path(__FILE__).parent_path();
namespace aniso_diff {

class ParameterReader : public Subscriptor {
public:
  std::string lbl = "unknown";
  ParameterReader(ParameterHandler &paramhandler) : prm(paramhandler){};
  void read_parameters(std::string param_fpath) {
    // Check file exists
    if (!fs::exists(param_fpath)) {
      fs::path qualified_path = this_dir / fs::path(param_fpath);
      Assert(fs::exists(qualified_path), ExcIO("Couldn't find param file"));
      param_fpath = qualified_path.string();
    }
    this->lbl = fs::path(param_fpath).stem();
    declare_parameters();
    prm.parse_input(param_fpath);
  };

private:
  void declare_parameters() {
    prm.enter_subsection("General");
    { prm.declare_entry("label", this->lbl); }
    prm.leave_subsection();
    prm.enter_subsection("FEM");
    {
      prm.declare_entry("degree", "2", Patterns::Integer(0),
                        "Basis function degree");
      prm.declare_entry("nref", "3", Patterns::Integer(1),
                        "Number of times to refine square mesh");
    }
    prm.leave_subsection();

    prm.enter_subsection("solver");
    {
      prm.declare_entry("type", "Bicgstab");
      prm.declare_entry("preconditioner", "jacobi");
      prm.declare_entry("nit", "1000", Patterns::Integer(11),
                        "Max number of solver iterations");
      prm.declare_entry("rtol", "1e-10", Patterns::Double(0, 1),
                        "Relative solver tolerance");
    }
    prm.leave_subsection();

    prm.enter_subsection("Phys");
    {
      prm.declare_entry("alpha", "0.0", Patterns::Double(0),
                        "Magnetic field parameter");
      prm.declare_entry("m", "0.0", Patterns::Double(0),
                        "Magnetic field parameter");
      prm.declare_entry("epsilon", "0.0", Patterns::Double(0), "Dpar/Dperp");
    }
    prm.leave_subsection();
  };
  ParameterHandler &prm;
};

template <int dim> class AnisotropicDiffusion {
public:
  AnisotropicDiffusion(ParameterHandler &param, const int degree);
  ~AnisotropicDiffusion();

  void run();

private:
  void setup_geometry();
  void setup_dof_handler();
  void setup_system();
  void assemble_system();
  void solve_system();
  void output_results(const std::string &suffix = "",
                      const std::string &fbase = "aniso-diff") const;

  Triangulation<dim> triangulation;
  DoFHandler<dim, dim> dof_handler;
  SparsityPattern sparsity_pattern;
  FE_Q<dim> fe;
  SparseMatrix<double> system_matrix;
  Vector<double> solution;
  Vector<double> system_rhs;

  std::string run_label;
  int nref;
  int solver_nit;
  double solver_rtol;
  std::string preconditioner_type = "jacobi";
  std::string solver_type;

  double alpha;
  double eps;
  double m;

  // Boundary values function
  class LowX_Bdy_Func : public Function<dim> {
  public:
    virtual double
    value([[maybe_unused]] const Point<dim> &p,
          [[maybe_unused]] const unsigned int component = 0) const override {
      return 0.0;
    }
  };
  LowX_Bdy_Func lowx_bdy_func;

  class ExactSln : public Function<dim> {
  public:
    ExactSln(const double alpha_in, const double m_in, const double eps_in)
        : alpha(alpha_in), eps(eps_in), m(m_in){};
    const double alpha, eps, m;
    virtual double
    value([[maybe_unused]] const Point<dim> &p,
          [[maybe_unused]] const unsigned int component = 0) const override {
      return sin(M_PI * p[1] + this->alpha * (p[1] * p[1] - p[1]) *
                                   cos(this->m * M_PI * p[0])) +
             this->eps * cos(2 * M_PI * p[0]) * sin(M_PI * p[1]);
    }
  };

  class DiffTensor : public TensorFunction<2, dim> {
  public:
    DiffTensor(const double alpha_in, const double m_in, const double eps_in)
        : TensorFunction<2, dim>(), alpha(alpha_in), kdiff(1.0 - eps_in),
          kperp(eps_in), m(m_in){};
    const double alpha, kdiff, kperp, m;
    virtual Tensor<2, dim> value(const Point<dim> &p) const override {
      Tensor<2, dim> out;
      const double x = p[0], y = p[1];
      const double bhat_x = alpha * (2 * y - 1) * cos(m * M_PI * x) + M_PI;
      const double bhat_y = M_PI * alpha * m * (y * y - y) * sin(m * M_PI * x);
      const double norm = std::sqrt(bhat_x * bhat_x + bhat_y * bhat_y),
                   bx = bhat_x / norm, by = bhat_y / norm;
      out[0][0] = kdiff * bx * bx + kperp;
      out[0][1] = kdiff * bx * by;
      out[1][0] = kdiff * bx * by;
      out[1][1] = kdiff * by * by + kperp;
      return out;
    }
  };

  // Right-hand side function
  struct RHSFunction : public Function<dim> {
    RHSFunction(const double alpha_in, const double m_in, const double eps_in)
        : alpha(alpha_in), epsilon(eps_in), m(m_in){};

  public:
    const double alpha, epsilon, m;
    virtual double
    value([[maybe_unused]] const Point<dim> &p,
          [[maybe_unused]] const unsigned int component = 0) const override {
      const double x = p[0];
      const double y = p[1];
      // Add source term func here
      return -std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                 (1 - epsilon) * (2 * y - 1) * (std::pow(y, 2) - y) *
                 (M_PI * epsilon * std::cos(2 * M_PI * x) * std::cos(M_PI * y) +
                  (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y)) *
                 std::pow(std::sin(M_PI * m * x), 2) /
                 (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      std::pow(std::pow(y, 2) - y, 2) *
                      std::pow(std::sin(M_PI * m * x), 2) +
                  std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                           2)) +
             2 * M_PI * std::pow(alpha, 2) * m * (1 - epsilon) *
                 (std::pow(y, 2) - y) *
                 (-M_PI * alpha * m * (std::pow(y, 2) - y) *
                      std::sin(M_PI * m * x) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) -
                  2 * M_PI * epsilon * std::sin(2 * M_PI * x) *
                      std::sin(M_PI * y)) *
                 std::sin(M_PI * m * x) * std::cos(M_PI * m * x) /
                 (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      std::pow(std::pow(y, 2) - y, 2) *
                      std::pow(std::sin(M_PI * m * x), 2) +
                  std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                           2)) +
             std::pow(M_PI, 2) * alpha * std::pow(m, 2) * (1 - epsilon) *
                 (std::pow(y, 2) - y) *
                 (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                 (M_PI * epsilon * std::cos(2 * M_PI * x) * std::cos(M_PI * y) +
                  (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y)) *
                 std::cos(M_PI * m * x) /
                 (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      std::pow(std::pow(y, 2) - y, 2) *
                      std::pow(std::sin(M_PI * m * x), 2) +
                  std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                           2)) +
             M_PI * alpha * m * (1 - epsilon) * (2 * y - 1) *
                 (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                 (-M_PI * alpha * m * (std::pow(y, 2) - y) *
                      std::sin(M_PI * m * x) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) -
                  2 * M_PI * epsilon * std::sin(2 * M_PI * x) *
                      std::sin(M_PI * y)) *
                 std::sin(M_PI * m * x) /
                 (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      std::pow(std::pow(y, 2) - y, 2) *
                      std::pow(std::sin(M_PI * m * x), 2) +
                  std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                           2)) +
             M_PI * alpha * m * (1 - epsilon) * (std::pow(y, 2) - y) *
                 (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                 (M_PI * epsilon * std::cos(2 * M_PI * x) * std::cos(M_PI * y) +
                  (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y)) *
                 (-2 * std::pow(M_PI, 3) * std::pow(alpha, 2) * std::pow(m, 3) *
                      std::pow(std::pow(y, 2) - y, 2) * std::sin(M_PI * m * x) *
                      std::cos(M_PI * m * x) +
                  2 * M_PI * alpha * m * (2 * y - 1) *
                      (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                      std::sin(M_PI * m * x)) *
                 std::sin(M_PI * m * x) /
                 std::pow(
                     std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                             std::pow(std::pow(y, 2) - y, 2) *
                             std::pow(std::sin(M_PI * m * x), 2) +
                         std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) +
                                      M_PI,
                                  2),
                     2) +
             2 * M_PI * alpha * m * (1 - epsilon) * (std::pow(y, 2) - y) *
                 (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                 (-M_PI * alpha * m * (2 * y - 1) * std::sin(M_PI * m * x) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) +
                  M_PI * alpha * m * (std::pow(y, 2) - y) *
                      (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                      std::sin(M_PI * m * x) *
                      std::sin(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) -
                  2 * std::pow(M_PI, 2) * epsilon * std::sin(2 * M_PI * x) *
                      std::cos(M_PI * y)) *
                 std::sin(M_PI * m * x) /
                 (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      std::pow(std::pow(y, 2) - y, 2) *
                      std::pow(std::sin(M_PI * m * x), 2) +
                  std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                           2)) +
             M_PI * alpha * m * (1 - epsilon) * (std::pow(y, 2) - y) *
                 (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                 (-M_PI * alpha * m * (std::pow(y, 2) - y) *
                      std::sin(M_PI * m * x) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) -
                  2 * M_PI * epsilon * std::sin(2 * M_PI * x) *
                      std::sin(M_PI * y)) *
                 (-std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      (4 * y - 2) * (std::pow(y, 2) - y) *
                      std::pow(std::sin(M_PI * m * x), 2) -
                  4 * alpha *
                      (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                      std::cos(M_PI * m * x)) *
                 std::sin(M_PI * m * x) /
                 std::pow(
                     std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                             std::pow(std::pow(y, 2) - y, 2) *
                             std::pow(std::sin(M_PI * m * x), 2) +
                         std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) +
                                      M_PI,
                                  2),
                     2) +
             (epsilon +
              (1 - epsilon) *
                  std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                           2) /
                  (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                       std::pow(std::pow(y, 2) - y, 2) *
                       std::pow(std::sin(M_PI * m * x), 2) +
                   std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                            2))) *
                 (-std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      std::pow(std::pow(y, 2) - y, 2) *
                      std::pow(std::sin(M_PI * m * x), 2) *
                      std::sin(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) -
                  std::pow(M_PI, 2) * alpha * std::pow(m, 2) *
                      (std::pow(y, 2) - y) * std::cos(M_PI * m * x) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) -
                  4 * std::pow(M_PI, 2) * epsilon * std::sin(M_PI * y) *
                      std::cos(2 * M_PI * x)) +
             (M_PI * epsilon * std::cos(2 * M_PI * x) * std::cos(M_PI * y) +
              (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                  std::cos(alpha * (std::pow(y, 2) - y) *
                               std::cos(M_PI * m * x) +
                           M_PI * y)) *
                 (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      (1 - epsilon) * (4 * y - 2) * (std::pow(y, 2) - y) *
                      std::pow(std::sin(M_PI * m * x), 2) /
                      (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                           std::pow(std::pow(y, 2) - y, 2) *
                           std::pow(std::sin(M_PI * m * x), 2) +
                       std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) +
                                    M_PI,
                                2)) +
                  std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                      (1 - epsilon) * std::pow(std::pow(y, 2) - y, 2) *
                      (-std::pow(M_PI, 2) * std::pow(alpha, 2) *
                           std::pow(m, 2) * (4 * y - 2) * (std::pow(y, 2) - y) *
                           std::pow(std::sin(M_PI * m * x), 2) -
                       4 * alpha *
                           (alpha * (2 * y - 1) * std::cos(M_PI * m * x) +
                            M_PI) *
                           std::cos(M_PI * m * x)) *
                      std::pow(std::sin(M_PI * m * x), 2) /
                      std::pow(std::pow(M_PI, 2) * std::pow(alpha, 2) *
                                       std::pow(m, 2) *
                                       std::pow(std::pow(y, 2) - y, 2) *
                                       std::pow(std::sin(M_PI * m * x), 2) +
                                   std::pow(alpha * (2 * y - 1) *
                                                    std::cos(M_PI * m * x) +
                                                M_PI,
                                            2),
                               2)) +
             (-M_PI * alpha * m * (std::pow(y, 2) - y) *
                  std::sin(M_PI * m * x) *
                  std::cos(alpha * (std::pow(y, 2) - y) *
                               std::cos(M_PI * m * x) +
                           M_PI * y) -
              2 * M_PI * epsilon * std::sin(2 * M_PI * x) *
                  std::sin(M_PI * y)) *
                 (-2 * M_PI * alpha * m * (1 - epsilon) * (2 * y - 1) *
                      (alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI) *
                      std::sin(M_PI * m * x) /
                      (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                           std::pow(std::pow(y, 2) - y, 2) *
                           std::pow(std::sin(M_PI * m * x), 2) +
                       std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) +
                                    M_PI,
                                2)) +
                  (1 - epsilon) *
                      std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) +
                                   M_PI,
                               2) *
                      (-2 * std::pow(M_PI, 3) * std::pow(alpha, 2) *
                           std::pow(m, 3) * std::pow(std::pow(y, 2) - y, 2) *
                           std::sin(M_PI * m * x) * std::cos(M_PI * m * x) +
                       2 * M_PI * alpha * m * (2 * y - 1) *
                           (alpha * (2 * y - 1) * std::cos(M_PI * m * x) +
                            M_PI) *
                           std::sin(M_PI * m * x)) /
                      std::pow(std::pow(M_PI, 2) * std::pow(alpha, 2) *
                                       std::pow(m, 2) *
                                       std::pow(std::pow(y, 2) - y, 2) *
                                       std::pow(std::sin(M_PI * m * x), 2) +
                                   std::pow(alpha * (2 * y - 1) *
                                                    std::cos(M_PI * m * x) +
                                                M_PI,
                                            2),
                               2)) +
             (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                  (1 - epsilon) * std::pow(std::pow(y, 2) - y, 2) *
                  std::pow(std::sin(M_PI * m * x), 2) /
                  (std::pow(M_PI, 2) * std::pow(alpha, 2) * std::pow(m, 2) *
                       std::pow(std::pow(y, 2) - y, 2) *
                       std::pow(std::sin(M_PI * m * x), 2) +
                   std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                            2)) +
              epsilon) *
                 (2 * alpha * std::cos(M_PI * m * x) *
                      std::cos(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y) -
                  std::pow(M_PI, 2) * epsilon * std::sin(M_PI * y) *
                      std::cos(2 * M_PI * x) -
                  std::pow(alpha * (2 * y - 1) * std::cos(M_PI * m * x) + M_PI,
                           2) *
                      std::sin(alpha * (std::pow(y, 2) - y) *
                                   std::cos(M_PI * m * x) +
                               M_PI * y));
    }
  };
};

template <int dim>
AnisotropicDiffusion<dim>::AnisotropicDiffusion(ParameterHandler &param,
                                                const int degree)
    : dof_handler(triangulation), fe(degree) {

  param.enter_subsection("FEM");
  this->nref = param.get_integer("nref");
  param.leave_subsection();

  param.enter_subsection("solver");
  this->solver_nit = param.get_double("nit");
  this->solver_rtol = param.get_double("rtol");
  this->solver_type = param.get("type");
  this->preconditioner_type = param.get("preconditioner");
  std::transform(this->solver_type.begin(), this->solver_type.end(),
                 this->solver_type.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  this->preconditioner_type = param.get("preconditioner");
  param.leave_subsection();

  // Model params
  param.enter_subsection("Phys");
  this->alpha = param.get_double("alpha");
  this->m = param.get_double("m");
  this->eps = param.get_double("epsilon");
  param.leave_subsection();

  param.enter_subsection("General");
  this->run_label = param.get("label");
  param.leave_subsection();
}

template <int dim> AnisotropicDiffusion<dim>::~AnisotropicDiffusion() {}

template <int dim> void AnisotropicDiffusion<dim>::setup_geometry() {
  const double dim_min = 0.0, dim_max = 1.0;
  GridGenerator::hyper_cube(triangulation, dim_min, dim_max, true);
  triangulation.refine_global(this->nref);
}

template <int dim> void AnisotropicDiffusion<dim>::assemble_system() {
  const QGauss<dim> quadrature(fe.degree + 2);
  FEValues<dim> fe_values(fe, quadrature,
                          update_values | update_gradients | update_JxW_values |
                              update_quadrature_points);
  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  // Matrices to hold cell values
  FullMatrix<double> cell_matrix;
  cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const RHSFunction rhs_func(this->alpha, this->m, this->eps);
  Vector<double> cell_rhs(dof_handler.n_dofs());
  Vector<double> rhs_values(dof_handler.n_dofs());
  Tensor<2, dim> diff_tensor;
  const DiffTensor diff_tensor_func(this->alpha, this->m, this->eps);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    fe_values.reinit(cell);
    cell_matrix = 0;
    cell_rhs = 0;
    const auto JxW = fe_values.get_JxW_values();
    for (unsigned int q = 0; q < quadrature.size(); ++q) {
      const auto p = fe_values.quadrature_point(q);
      diff_tensor = diff_tensor_func.value(p);
      double rhs_value = rhs_func.value(p);
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) {
        auto gradi_u = fe_values.shape_grad(i, q);
        for (unsigned int j = 0; j < fe.dofs_per_cell; ++j) {
          auto gradj_u = fe_values.shape_grad(j, q);
          cell_matrix.add(i, j,
                          JxW[q] *
                              (diff_tensor[0][0] * gradi_u[0] * gradj_u[0] +
                               diff_tensor[0][1] * gradi_u[0] * gradj_u[1] +
                               diff_tensor[1][0] * gradi_u[1] * gradj_u[0] +
                               diff_tensor[1][1] * gradi_u[1] * gradj_u[1]));
        }
        // Right-hand side
        cell_rhs[i] += fe_values.shape_value(i, q) * rhs_value * JxW[q];
      }
    }

    std::vector<types::global_dof_index> local_dof_indices(
        fe.n_dofs_per_cell());
    cell->get_dof_indices(local_dof_indices);
    // Write cell values into system matrices
    for (auto ii : fe_values.dof_indices()) {
      for (auto jj : fe_values.dof_indices()) {
        system_matrix.add(local_dof_indices[ii], local_dof_indices[jj],
                          cell_matrix(ii, jj));
      }
      system_rhs(local_dof_indices[ii]) = cell_rhs(ii);
    }
  }

  std::map<types::global_dof_index, double> boundary_values;
  constexpr int left_bid = 0, right_bid = 1, bottom_bid = 2, top_bid = 3;
  // top, bottom boundaries homogeneous Dirichlet
  for (auto &bid : {top_bid, bottom_bid}) {
    VectorTools::interpolate_boundary_values(
        dof_handler, bid, ZeroFunction<dim>(), boundary_values);
  }
  for (auto &bid : {left_bid, right_bid}) {
    VectorTools::interpolate_boundary_values(
        dof_handler, bid, ExactSln(this->alpha, this->m, this->eps),
        boundary_values);
  }
  // Top, bottom boundaries are not set explicitly => homgeneous Neumann
  MatrixTools::apply_boundary_values(boundary_values, this->system_matrix,
                                     this->solution, this->system_rhs);
}

template <int dim> void AnisotropicDiffusion<dim>::setup_system() {
  this->dof_handler.distribute_dofs(this->fe);
  auto ndofs = this->dof_handler.n_dofs();

  DynamicSparsityPattern dsp(ndofs);
  DoFTools::make_sparsity_pattern(this->dof_handler, dsp);
  this->sparsity_pattern.copy_from(dsp);
  this->system_matrix.reinit(this->sparsity_pattern);
  this->solution.reinit(ndofs);
  this->system_rhs.reinit(ndofs);
}

template <int dim> void AnisotropicDiffusion<dim>::solve_system() {
  SolverControl solver_control(this->solver_nit,
                               this->solver_rtol * this->system_rhs.l2_norm());

  SolverSelector<Vector<double>> solver(this->solver_type, solver_control);

  PreconditionSelector<SparseMatrix<double>, Vector<double>> preconditioner(
      this->preconditioner_type, 1.);
  preconditioner.use_matrix(this->system_matrix);

  solver.solve(this->system_matrix, this->solution, this->system_rhs,
               preconditioner);
}

template <int dim>
void AnisotropicDiffusion<dim>::output_results(const std::string &suffix,
                                               const std::string &fbase) const {

  Vector<double> exact_solution;
  exact_solution.reinit(this->dof_handler.n_dofs());
  VectorTools::interpolate(
      dof_handler, ExactSln(this->alpha, this->m, this->eps), exact_solution);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(this->dof_handler);
  data_out.add_data_vector(this->solution, "solution");
  data_out.add_data_vector(exact_solution, "exact");

  data_out.build_patches();

  fs::path fpath = this_dir / fs::path("output") /
                   fs::path(this->run_label + "_" + fbase + suffix + ".vtk");
  std::ofstream output(fpath);
  data_out.write_vtk(output);
  deallog << " Wrote output to" << fpath;
}

template <int dim> void AnisotropicDiffusion<dim>::run() {
  setup_geometry();
  setup_system();
  assemble_system();
  solve_system();
  output_results();
}
} // namespace aniso_diff

int main(int argc, char *argv[]) {
  deallog.depth_console(2);

  constexpr int nargs_expected = 2;
  if (argc != 2) {
    std::cout << "Expected " << nargs_expected << " args, but got " << argc
              << std::endl;
    return 1;
  }

  ParameterHandler param_handler;
  aniso_diff::ParameterReader param_reader(param_handler);
  param_reader.read_parameters(argv[1]);

  // Set basis degree from params
  param_handler.enter_subsection("FEM");
  int degree = param_handler.get_integer("degree");
  param_handler.leave_subsection();
  aniso_diff::AnisotropicDiffusion<2> prob(param_handler, degree);
  prob.run();

  return 0;
}
