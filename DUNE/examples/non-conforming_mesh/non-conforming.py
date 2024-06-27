from dune.alugrid import aluGrid, aluConformGrid, aluCubeGrid, aluSimplexGrid
from dune.grid import reader
import meshio
import os.path
import sys

# Some useful locations
repo_root = os.path.abspath(os.path.dirname(__file__) + "../../../../")
this_dir = os.path.abspath(os.path.dirname(__file__))


def convert(
    input_fname,
    grid_from_mesh_data=False,
    read_args={},
    write_dgf_and_read_back=False,
    write_msh=False,
    write_vtk=True,
):
    mesh = read_mesh(input_fname, **read_args)
    fname_base = os.path.splitext(input_fname)[0]
    if write_msh:
        msh_path = write_to_msh(mesh, fname_base + "_out")
    if write_vtk:
        vtk_path = write_to_vtk(mesh, fname_base + "_out")
    if write_dgf_and_read_back:
        dgf_path = write_to_dgf(mesh, fname_base)
        grid = grid_from_dgf(dgf_path, aluCubeGrid)
        # grid = grid_from_dgf(dgf_path,aluSimplexGrid)
        grid.plot(figsize=(5, 5))
    elif grid_from_mesh_data:
        grid = grid_from_mesh(mesh)
        grid.plot(figsize=(5, 5))


def read_mesh(fname, dir=this_dir, **read_args):
    fpath = os.path.join(dir, fname)
    print(f"Reading mesh from {fpath}")
    return meshio.read(fpath, **read_args)


def grid_from_mesh(mesh):
    # Construct grid from gmsh data; fails
    points, cells = mesh.points, mesh.cells_dict
    return aluConformGrid({"vertices": points, "simplices": cells})


def grid_from_dgf(fpath, grid_func):
    domain = (reader.dgf, fpath)
    return grid_func(domain, dimgrid=2)


def write_to_dgf(mesh, basename, dir=this_dir):
    # Add scripts dir to path and import conversion util
    sys.path.append(os.path.join(repo_root, "DUNE/scripts/"))
    from gmsh2DGF import gmsh2DGF

    points, cells = mesh.points, mesh.cells_dict
    dgf = gmsh2DGF(points, cells)
    outpath = os.path.join(dir, f"{basename}.dgf")
    with open(outpath, "w") as file:
        file.write(dgf)
        print(f"Wrote DGF to {outpath}")
    return outpath


def write_to_msh(mesh, basename, dir=this_dir):
    outpath = os.path.join(dir, f"{basename}.msh")
    meshio.write(outpath, mesh, file_format="gmsh", binary=False)
    print(f"Wrote MSH to {outpath}")


def write_to_vtk(mesh, basename, dir=this_dir):
    outpath = os.path.join(dir, f"{basename}.vtk")
    meshio.write(outpath, mesh, binary=False)
    print(f"Wrote VTK to {outpath}")


convert(
    "domain_4_byhand.msh",
    grid_from_mesh_data=False,
    write_dgf_and_read_back=True,
    write_msh=True,
    write_vtk=False,
    read_args=dict(file_format="gmsh"),
)
