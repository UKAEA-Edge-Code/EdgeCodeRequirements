import meshio
import os.path
import sys

# Add scripts dir to path and import conversion util
repo_root = os.path.abspath(os.path.dirname(__file__) + "../../../../../")
sys.path.append(os.path.join(repo_root, "DUNE/scripts/"))
from gmsh2DGF import gmsh2DGF

# Read mesh file from Firedrake eg and extract points, cells
inpath = os.path.join(
    repo_root,
    "Firedrake/examples/scripts/aniso_diffusion/aniso_diffusion_DeluzetNarski.msh",
)
mesh = meshio.read(inpath)
points, cells = mesh.points, mesh.cells_dict

# Generate dgf and output to this dir
outpath = os.path.join(os.path.abspath(os.path.dirname(__file__)), "unit_square.dgf")
dgf = gmsh2DGF(points, cells)
with open(outpath, "w") as file:
    file.write(dgf)
    print(f"Wrote mesh to {outpath}")
