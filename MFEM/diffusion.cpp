// Adapated from MFEM Example 1
//
// Compile with: make
//

#include "mfem.hpp"
#include "admfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;
using namespace mfem;

// This should be defined in MFEM, but doesn't seem to be in the version I have installed
typedef double real_t;


template <typename T> void diffusivity_coefficients(const Vector &params, const T &input, T &output) {
  const auto alpha = params[0], m = params[1], epsilon = params[2];
  const auto x = input[0], y = input[1];
  const auto kpar = 1.0, kperp = epsilon, kdiff = kpar - kperp;
  const auto x_comp = alpha * (2 * y - 1) * cos(m * M_PI * x) + M_PI,
    y_comp = M_PI * alpha * m * (y * y - y) * sin(m * M_PI * x);
  const auto denom = sqrt(x_comp * x_comp + y_comp * y_comp), bx = x_comp / denom, by = y_comp / denom;
  output[0] = kdiff * bx * bx + kperp;
  output[1] = kdiff * bx * by;
  output[2] = kdiff * bx * by;
  output[3] = kdiff * by * by + kperp;
}

template <typename ReturnType, typename ParamVector, typename StateVector,
          int state_size, int param_size>
class SolutionFunctor {
public:
  ReturnType operator()(const ParamVector &params, const StateVector &input) const {
    const auto alpha = params[0], m = params[1], epsilon = params[2];
    const auto x = input[0], y = input[1];
    return sin(M_PI * y + alpha * (y * y - y) * cos(m * M_PI * x)) + epsilon * cos(2 * M_PI * x) * sin(M_PI * y);
  }
};

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "square.msh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   bool fa = false;
   real_t epsilon = 1e-1;
   real_t alpha = 2.;
   real_t m = 1.;
   const char *device_config = "cpu";
   bool visualization = true;
   bool paraview = false;
   bool algebraic_ceed = false;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&fa, "-fa", "--full-assembly", "-no-fa",
                  "--no-full-assembly", "Enable Full Assembly.");
   args.AddOption(&epsilon, "-eps", "--epsilon",
                  "Ratio between perpendicular and parallel diffusivity."
                  );
   args.AddOption(&alpha, "-a", "--alpha", "Alpha parameter for magnetic field.");
   args.AddOption(&m, "-mp", "--m-param", "m parameter for magnetic field.");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
#ifdef MFEM_USE_CEED
   args.AddOption(&algebraic_ceed, "-a", "--algebraic", "-no-a", "--no-algebraic",
                  "Use algebraic Ceed solver");
#endif
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&paraview, "-paraview", "--paraview-datafiles", "-no-paraview",
                  "--no-paraview-datafiles",
                  "Save data files for ParaView visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   Vector field_params({alpha, m, epsilon});

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 4. Set the desired solution and manufacture a forcing term that will
   //    produce it.
   const SolutionFunctor<real_t, Vector, Vector, 2, 3> real_solution;
   FunctionCoefficient expected([& field_params, & real_solution](const Vector & x) -> real_t {
     return real_solution(field_params, x);
   }
     );
   const SolutionFunctor<ad::ADFloatType, Vector, ad::ADVectorType, 2, 3> differentiable_solution;
   // Use automatic differentiation to compute the forcing term, see https://mfem.org/autodiff/
   VectorFuncAutoDiff<1, 2, 3> solution_grad([& differentiable_solution](const Vector &params, const ad::ADVectorType &input, ad::ADVectorType &output) {
     output[0] = differentiable_solution(params, input);
   });
   QFunctionAutoDiff<SolutionFunctor, 2, 3> solution_hessian;
   VectorFuncAutoDiff<4, 2, 3> diffusivity_grad(diffusivity_coefficients<ad::ADVectorType>);
   FunctionCoefficient forcing(
                               [& field_params, & solution_grad, & solution_hessian,
        & diffusivity_grad](const Vector &x) -> real_t {
                                  // The function calls below don't
                                  // take const arguments, although
                                  // pretty sure they could have been
                                  // written to
                                 Vector x_var(x), diff(4);
                                 DenseMatrix sol_grad(1, 2), sol_hess(2, 2), diff_jac(4, 2);
                                 diffusivity_coefficients(field_params, x, diff);
                                 solution_grad.Jacobian(field_params, x_var, sol_grad);
                                 solution_hessian.Hessian(field_params, x_var, sol_hess);
                                 diffusivity_grad.Jacobian(field_params, x_var, diff_jac);
                                 real_t dj00 = diff_jac(0, 0), dj10 = diff_jac(1, 0), dj21 = diff_jac(2, 1), dj31 = diff_jac(3,1), sg00 = sol_grad(0,0), sg01 = sol_grad(0, 1), sh00 = sol_hess(0, 0), sh01 = sol_hess(0, 1), sh10 = sol_hess(1, 0), sh11 = sol_hess(1, 1);
                                 real_t result= -(
                                          diff_jac(0, 0) * sol_grad(0, 0) + diff[0] * sol_hess(0, 0) +
                                          diff_jac(1, 0) * sol_grad(0, 1) + diff[1] * sol_hess(0, 1) +
                                          diff_jac(2, 1) * sol_grad(0, 0) + diff[2] * sol_hess(1, 0) +
                                          diff_jac(3, 1) * sol_grad(0, 1) + diff[3] * sol_hess(1, 1));
                                 return result;
                               });

   // 5. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   bool delete_fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
      delete_fec = true;
   }
   else if (mesh.GetNodes())
   {
      fec = mesh.GetNodes()->OwnFEC();
      delete_fec = false;
      cout << "Using isoparametric FEs: " << fec->Name() << endl;
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
      delete_fec = true;
   }
   FiniteElementSpace fespace(&mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace.GetTrueVSize() << endl;

   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   Array<int> ess_tdof_list;
   Array<int> ess_bdr(mesh.bdr_attributes.Max());
   if (mesh.bdr_attributes.Size())
   {
      // Boundaries 11 and 13 from square.geo are Dirichlet, but MFEM
      // uses 0-indexing so subtract 1.
     // Degenerate? Try imposing Dirichlet on all boundaries
      ess_bdr = 1;
      ess_bdr[10] = 1;
      ess_bdr[12] = 1;
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
      // for (auto &&i : ess_tdof_list) {
      //   std::cout << i << std::endl;
      // }
   }

   // 7. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   LinearForm b(&fespace);
   b.AddDomainIntegrator(new DomainLFIntegrator(forcing));
   b.Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace);   
   x = 0.0;
   x.ProjectBdrCoefficient(expected, ess_bdr);

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   MatrixFunctionCoefficient diffusivity(2, [field_params](const Vector & x, DenseMatrix & mat){
     Vector coeffs(4);
     diffusivity_coefficients(field_params, x, coeffs);
     mat(0, 0) = coeffs[0];
     mat(0, 1) = coeffs[1];
     mat(1, 0) = coeffs[2];
     mat(1, 1) = coeffs[3];
   });

   BilinearForm a(&fespace);
   if (pa) { a.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   if (fa)
   {
      a.SetAssemblyLevel(AssemblyLevel::FULL);
      // Sort the matrix column indices when running on GPU or with OpenMP (i.e.
      // when Device::IsEnabled() returns true). This makes the results
      // bit-for-bit deterministic at the cost of somewhat longer run time.
      a.EnableSparseMatrixSorting(Device::IsEnabled());
   }
   a.AddDomainIntegrator(new DiffusionIntegrator(diffusivity));

   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   if (static_cond) { a.EnableStaticCondensation(); }
   a.Assemble();

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   cout << "Size of linear system: " << A->Height() << endl;

   // 11. Solve the linear system A X = B.
   if (!pa)
   {
#ifndef MFEM_USE_SUITESPARSE
      // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
      GSSmoother M((SparseMatrix&)(*A));
      PCG(*A, M, B, X, 1, 1000, 1e-12, 0.0);
#else
      // If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
      UMFPackSolver umf_solver;
      umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      umf_solver.SetOperator(*A);
      umf_solver.Mult(B, X);
#endif
   }
   else
   {
      if (UsesTensorBasis(fespace))
      {
         if (algebraic_ceed)
         {
            ceed::AlgebraicSolver M(a, ess_tdof_list);
            PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
         }
         else
         {
            OperatorJacobiSmoother M(a, ess_tdof_list);
            PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
         }
      }
      else
      {
         CG(*A, B, X, 1, 400, 1e-12, 0.0);
      }
   }

   // 12. Recover the solution as a finite element grid function.
   a.RecoverFEMSolution(X, b, x);
   const real_t l2 = x.ComputeL2Error(expected);
   cout << "L2 error norm: " << l2 << endl;

   // 13. Save the mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m diffusion.mesh -g sol.gf".
   ofstream mesh_ofs("diffusion.mesh");
   mesh_ofs.precision(8);
   mesh.Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
   sol_ofs.precision(8);
   x.Save(sol_ofs);

   // 14. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << x << flush;
   }

   // 14. Save the solution in paraview format.
   if (paraview) {
     ParaViewDataCollection pd("diffusion", &mesh);
     pd.SetPrefixPath("ParaView");
     pd.RegisterField("solution", &x);
     GridFunction expected_x(&fespace);
     expected_x.ProjectCoefficient(expected);
     pd.RegisterField("expected", &expected_x);
     pd.SetLevelsOfDetail(order);
     pd.SetDataFormat(VTKFormat::BINARY);
     pd.SetHighOrderOutput(true);
     pd.SetCycle(0);
     pd.SetTime(0.);
     pd.Save();
   }

   // 15. Free the used memory.
   if (delete_fec)
   {
      delete fec;
   }

   return 0;
}
