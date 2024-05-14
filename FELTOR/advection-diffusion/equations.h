// equations.h

#pragma once

#include "dg/algorithm.h"

namespace myproject
{

    template <class Geometry, class Matrix, class Container>
    struct Equations
    {
        Equations(const Geometry &grid, dg::file::WrappedJsonValue &js) : m_phi(dg::evaluate(dg::zero, grid)),
                                                                          m_old_phi(2, m_phi),
                                                                          m_adv(grid), // use grid's boundary conditions
                                                                          m_lapM(grid)
        {
            m_v = {m_phi, m_phi};
            m_pcg.construct(m_phi, grid.size());
            m_eps_pol = js["elliptic"].get("eps_pol", 1e-6).asDouble();
            m_nu = js["physical"].get("nu", 1e-7).asDouble();
            m_centered[0] = dg::create::dx(grid, grid.bcx(), dg::centered);
            m_centered[1] = dg::create::dy(grid, grid.bcy(), dg::centered);
        }
        // accessors for diag.h
        const Container &potential() const { return m_phi; }
        // We implement advection diffusion equations:
        void operator()(double t, const Container &omega, Container &omegaDot)
        {

            // Solve potential equation
            m_old_phi.extrapolate(t, m_phi);
            // For demonstration we here use the simple unpreconditioned PCG
            // solver. In real code it is highly recommended to use nested_iterations instead
            // See "Solvers" chapter of this guide
            m_pcg.solve(m_lapM, m_phi, omega, 1., m_lapM.weights(), m_eps_pol);
            m_old_phi.update(t, m_phi);

            // add advection term
            dg::blas2::symv(-1., m_centered[1], m_phi, 0., m_v[0]);
            dg::blas2::symv(+1., m_centered[0], m_phi, 0., m_v[1]);
            dg::blas1::plus(m_v[0], 1.);
            m_adv.upwind(-1., m_v[0], m_v[1], omega, 0., omegaDot);

            // add diffusion
            dg::blas2::symv(-m_nu, m_lapM, omega, 1., omegaDot);
        }

    private:
        Container m_phi;
        dg::Extrapolation<Container> m_old_phi;
        dg::Advection<Geometry, Matrix, Container> m_adv;
        dg::Elliptic<Geometry, Matrix, Container> m_lapM;
        dg::PCG<Container> m_pcg;
        std::array<Matrix, 2> m_centered;
        std::array<Container, 2> m_v;
        double m_eps_pol;
        double m_nu;
    };

} // namespace myproject