// diag.h
#pragma once
#include "equations.h"

namespace myproject
{

    struct Variables
    {
        Equations<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec> &rhs;
        const dg::x::CartesianGrid2d &grid;
        const dg::x::DVec &omega;
    };

    struct Record
    {
        std::string name;      // variable name in the output file
        std::string long_name; // longer description as an attribute
        std::function<void(dg::x::DVec &, Variables &)> function;
        // function that generates the data points for the variable
    };

    // time - independent output (only called once)
    std::vector<Record> diagnostics2d_static_list = {
        {"xc", "x-coordinate in Cartesian coordinate system",
         [](dg::x::DVec &result, Variables &v)
         {
             result = dg::evaluate(dg::cooX2d, v.grid);
         }},
        {"yc", "y-coordinate in Cartesian coordinate system",
         [](dg::x::DVec &result, Variables &v)
         {
             result = dg::evaluate(dg::cooY2d, v.grid);
         }},
        {"weights", "Gaussian integration weights",
         [](dg::x::DVec &result, Variables &v)
         {
             result = dg::create::weights(v.grid);
         }},
        // ... extend here
    };

    // time - dependent output (called periodically)
    std::vector<Record> diagnostics2d_list = {
        {"vorticity", "Vorticity in 2d",
         [](dg::x::DVec &result, Variables &v)
         {
             dg::blas1::copy(v.omega, result);
         }},
        {"potential", "stream function",
         [](dg::x::DVec &result, Variables &v)
         {
             dg::blas1::copy(v.rhs.potential(), result);
         }},
        {"enstrophy", "Squared vorticity",
         [](dg::x::DVec &result, Variables &v)
         {
             // more complicated algorithms can be written here
             dg::blas1::pointwiseDot(v.omega, v.omega, result);
             dg::blas1::scal(result, 1. / 2.);
         }}
        // ... extend here
    };

    // ... write more lists, for example 1d diagnostics etc
    // each list requires a loop in the main program

} // namespace myproject