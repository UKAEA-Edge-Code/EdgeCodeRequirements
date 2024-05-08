// init.h
#pragma once

namespace myproject
{

dg::x::DVec initial_conditions(
    const dg::x::CartesianGrid2d& grid,
    dg::file::WrappedJsonValue init)
{
    dg::x::HVec omega;
    std::string initial = init.get("type", "lamb").asString();
    if( initial == "zero")
    {
        omega = dg::evaluate( dg::zero, grid);
    }
    else if( initial == "lamb")
    {
        double posX = init.get("posX", 0.5).asDouble();
        double posY = init.get("posY", 0.8).asDouble();
        double R =    init.get("sigma", 0.1).asDouble();
        double U =    init.get("velocity", 1).asDouble();
        dg::Lamb lamb( posX*grid.lx(), posY*grid.ly(), R, U);
        omega = dg::evaluate ( lamb, grid);
    }
    else if (initial == "gaussian")
    {
        double posX = init.get("posX", 0.2).asDouble();
        double posY = init.get("posY", 0.5).asDouble();
        double R =    init.get("sigma", 1).asDouble();
        dg::Gaussian gaussian( posX*grid.lx(), posY*grid.ly(), R, R, 1);
        omega = dg::evaluate ( gaussian, grid);
    }
    // ... implement more initial conditions here
    else
        throw dg::Error( dg::Message() << "Initial condition "
                        <<initial<<" not recognized!");
    return dg::construct<dg::x::DVec>(omega);
}

} //namespace myproject