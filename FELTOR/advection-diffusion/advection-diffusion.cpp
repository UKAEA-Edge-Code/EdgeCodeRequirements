// myprogram.cpp
#include <iostream>
#include <iomanip>

#ifdef WITH_MPI
#include <mpi.h>
#endif // WITH_MPI

#include "dg/algorithm.h"
#include "dg/file/file.h"

#include "equations.h"
#include "init.h"
#include "diag.h"

int main(int argc, char *argv[])
{
#ifdef WITH_MPI
    dg::mpi_init(argc, argv);
    MPI_Comm comm;
    dg::mpi_init2d(dg::DIR, dg::PER, comm, std::cin, true);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif // WITH_MPI
//
//
//Read in input
    dg::file::WrappedJsonValue js(dg::file::error::is_throw);
    if (argc != 3)
    {
        DG_RANK0 std::cerr << "ERROR: Wrong number of arguments!\nUsage: "
                           << argv[0] << " [input.json] [output.nc]\n \n"
                           << std::endl;
        dg::abort_program();
    }
    try
    {
        dg::file::file2Json(argv[1], js.asJson(),
                            dg::file::comments::are_discarded, dg::file::error::is_throw);
    }
    catch (std::exception &e)
    {
        DG_RANK0 std::cerr << "ERROR in input file " << argv[1] << std::endl;
        DG_RANK0 std::cerr << e.what() << std::endl;
        dg::abort_program();
    }
    //DG_RANK0 std::cout << js.asJson() << std::endl;
//
//
//
// Construct grid
    unsigned n, Nx, Ny;
    double x0, x1, y0, y1;
    dg::bc bcx, bcy;
    try
    {
        dg::file::WrappedJsonValue grid = js["grid"];
        n = grid.get("n", 3).asUInt();
        Nx = grid.get("Nx", 48).asUInt();
        Ny = grid.get("Ny", 48).asUInt();
        x0 = grid["x"].get(0u, 0.).asDouble();
        x1 = grid["x"].get(1u, 1.).asDouble();
        y0 = grid["y"].get(0u, 0.).asDouble();
        y1 = grid["y"].get(1u, 1.).asDouble();
        bcx = dg::str2bc(grid["bc"].get(0u, "DIR").asString());
        bcy = dg::str2bc(grid["bc"].get(1u, "PER").asString());
    }
    catch (std::exception &error)
    {
        DG_RANK0 std::cerr << "Error in input file " << argv[1] << std::endl;
        DG_RANK0 std::cerr << error.what() << std::endl;
        dg::abort_program();
    }
    dg::x::CartesianGrid2d grid(x0, x1, y0, y1, n, Nx, Ny, bcx, bcy
// The MPI version of CartesianGrid2d needs a communicator:
#ifdef WITH_MPI
                                ,
                                comm
#endif // WITH_MPI
    );
//
//
//
// Construct Equations
    myproject::Equations<dg::x::CartesianGrid2d, dg::x::DMatrix,
                         dg::x::DVec>
        rhs(grid, js);

//
// Construct initial condition
    dg::x::DVec omega;
    try
    {
        omega = myproject::initial_conditions(grid, js["init"]);
    }
    catch (std::exception &error)
    {
        DG_RANK0 std::cerr << "Error in input file " << argv[1] << std::endl;
        DG_RANK0 std::cerr << error.what() << std::endl;
        dg::abort_program();
    }

//
// Construct timestepper
    std::string tableau;
    double rtol, atol, time = 0.;
    try
    {
        rtol = js["timestepper"].get("rtol", 1e-5).asDouble();
        atol = js["timestepper"].get("atol", 1e-5).asDouble();
        tableau = js["timestepper"].get("tableau", "Bogacki-Shampine-4-2-3").asString();
    }
    catch (std::exception &error)
    {
        DG_RANK0 std::cerr << "Error in input file " << argv[1] << std::endl;
        DG_RANK0 std::cerr << error.what() << std::endl;
        dg::abort_program();
    }
    dg::Adaptive<dg::ERKStep<dg::x::DVec>> adapt(tableau, omega);
    dg::AdaptiveTimeloop<dg::x::DVec> timeloop(adapt, rhs,
                                               dg::pid_control, dg::l2norm, rtol, atol);
    myproject::Variables var = {rhs, grid, omega};
    // trigger first computation of potential
    {
        dg::x::DVec temp = omega;
        rhs(0., omega, temp);
    }

//
// Create netcdf file
    dg::file::NC_Error_Handle err;
    int ncid = -1;
    try
    {
        DG_RANK0 err = nc_create(argv[2], NC_NETCDF4 | NC_CLOBBER, &ncid);
    }
    catch (std::exception &e)
    {
        DG_RANK0 std::cerr << "ERROR creating file " << argv[1] << std::endl;
        DG_RANK0 std::cerr << e.what() << std::endl;
        dg::abort_program();
    }
    std::map<std::string, std::string> att;
    att["title"] = "Output file of myproject/myprogram.cpp";
    att["Conventions"] = "CF-1.8";
    /// Get local time and begin file history
    auto ttt = std::time(nullptr);

    std::ostringstream oss;
    /// time string  + program-name + args
    oss << std::put_time(std::localtime(&ttt), "%F %T %Z");
    for (int i = 0; i < argc; i++)
        oss << " " << argv[i];
    att["history"] = oss.str();
    att["comment"] = "Find more info in myproject/documentation.tex";
    att["source"] = "FELTOR";
    att["git-hash"] = GIT_HASH;
    att["git-branch"] = GIT_BRANCH;
    att["compile-time"] = COMPILE_TIME;
    att["references"] = "https://github.com/myname/myproject";
    // Here we put the inputfile as a string without comments so that it can be read later by another parser
    att["inputfile"] = js.asJson().toStyledString();
    for (auto pair : att)
        DG_RANK0 err = nc_put_att_text(ncid, NC_GLOBAL,
                                       pair.first.data(), pair.second.size(), pair.second.data());

//
// Set up Output                                       
    unsigned n_out = js["output"]["n"].asUInt(3);
    unsigned Nx_out = js["output"]["Nx"].asUInt(48);
    unsigned Ny_out = js["output"]["Ny"].asUInt(48);

    dg::x::CartesianGrid2d grid_out(x0, x1, y0, y1,
                                    n_out, Nx_out, Ny_out, bcx, bcy
#ifdef WITH_MPI
                                    ,
                                    comm
#endif // WITH_MPI
    );
    dg::x::IHMatrix projection = dg::create::interpolation(grid_out, grid);

//
//Define Dimensions and Variables
    int dim_ids[3], tvarID;
    // the dimensions are the ones of grid_out!
    err = dg::file::define_dimensions(ncid, dim_ids, &tvarID, grid_out,
                                      {"time", "y", "x"});

    std::map<std::string, int> id1d, id3d;
    for (auto &record : myproject::diagnostics2d_list)
    {
        std::string name = record.name;
        std::string long_name = record.long_name;
        id3d[name] = 0;
        DG_RANK0 err = nc_def_var(ncid, name.data(), NC_DOUBLE, 3, dim_ids,
                                  &id3d.at(name));
        DG_RANK0 err = nc_put_att_text(ncid, id3d.at(name), "long_name",
                                       long_name.size(), long_name.data());
        // and the 1d fields (our idea is to just volume integrate the 2d fields
        name = name + "_1d";
        long_name = long_name + " (Volume integrated)";
        id1d[name] = 0;
        DG_RANK0 err = nc_def_var(ncid, name.data(), NC_DOUBLE, 1, &dim_ids[0],
                                  &id1d.at(name));
        DG_RANK0 err = nc_put_att_text(ncid, id1d.at(name), "long_name",
                                       long_name.size(), long_name.data());
    }

//
//Output Static List
    dg::x::HVec resultH = dg::evaluate(dg::zero, grid);
    dg::x::HVec transferH = dg::evaluate(dg::zero, grid_out);
    dg::x::DVec resultD = transferH; // transfer to device
    for (auto &record : myproject::diagnostics2d_static_list)
    {
        std::string name = record.name;
        std::string long_name = record.long_name;
        int staticID = 0;
        DG_RANK0 err = nc_def_var(ncid, name.data(), NC_DOUBLE, 2, &dim_ids[1],
                                  &staticID);
        DG_RANK0 err = nc_put_att_text(ncid, staticID, "long_name",
                                       long_name.size(), long_name.data());
        record.function(resultD, var);
        dg::assign(resultD, resultH);
        dg::blas2::gemv(projection, resultH, transferH);
        dg::file::put_var_double(ncid, staticID, grid_out, transferH);
    }

//
//First File output
    dg::x::DVec volume = dg::create::volume(grid);
    size_t start = {0};
    size_t count = {1};
    for (auto &record : myproject::diagnostics2d_list)
    {
        record.function(resultD, var);
        double result = dg::blas1::dot(volume, resultD);
        dg::assign(resultD, resultH);
        dg::blas2::gemv(projection, resultH, transferH);
        // note that all processes call this function (for MPI)
        dg::file::put_vara_double(ncid, id3d.at(record.name), start,
                                  grid_out, transferH);
        // For the 1d output only the master thread needs to call the function
        DG_RANK0 err = nc_put_vara_double(ncid, id1d.at(record.name + "_1d"),
                                          &start, &count, &result);
    }
    DG_RANK0 err = nc_put_vara_double(ncid, tvarID, &start, &count, &time);
    DG_RANK0 err = nc_close(ncid);

//
//Timeloop
    double Tend = js["output"].get("tend", 1.0).asDouble();
    unsigned maxout = js["output"].get("maxout", 10).asUInt();
    double deltaT = Tend / (double)maxout;
    bool abort = false;
    for (unsigned u = 1; u <= maxout; u++)
    {

        try
        {
            // the documentation of dg::aTimeloop holds more details about how this construct works ...
            timeloop.integrate(time, omega, u * deltaT, omega,
                               u < maxout ? dg::to::at_least : dg::to::exact);

        }
        catch (std::exception &fail)
        {
            DG_RANK0 std::cerr << "ERROR in Timestepper\n";
            DG_RANK0 std::cerr << fail.what() << std::endl;
            DG_RANK0 std::cerr << "Writing last output and exit ..." << std::endl;
            abort = true;
        }
        start = u;
        DG_RANK0 err = nc_open(argv[2], NC_WRITE, &ncid);
        // First write the time variable
        DG_RANK0 err = nc_put_vara_double(ncid, tvarID, &start, &count, &time);
        for (auto &record : myproject::diagnostics2d_list)
        {
            record.function(resultD, var);
            double result = dg::blas1::dot(volume, resultD);
            dg::assign(resultD, resultH);
            dg::blas2::gemv(projection, resultH, transferH);
            dg::file::put_vara_double(ncid, id3d.at(record.name),
                                      start, grid_out, transferH);
            DG_RANK0 err = nc_put_vara_double(ncid, id1d.at(record.name + "_1d"),
                                              &start, &count, &result);
        }
        DG_RANK0 err = nc_close(ncid);
        if (abort)
            break;
    }

//
//Finalise
#ifdef WITH_MPI
    MPI_Finalize();
#endif // WITH_MPI
    return 0;
} // end of main