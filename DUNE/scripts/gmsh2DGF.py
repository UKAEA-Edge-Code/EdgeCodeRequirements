# %% [markdown]
#
# ## Converting gmsh to DGF
#
# A simple gmsh to DGF converter allowing to specify boundary conditions.
# A complete description of the DGF format is found [here](https://dune-project.org/doxygen/master/group__DuneGridFormatParser.html#details).
#
# %%
def gmsh2DGF(
    points,
    cells,
    bndDomain=None,
    bndSegments=None,
    periodic=None,
    dim=None,
    computeFirstIndex=True,
):
    """
    Parameter:
       points      array(list) of points of length dim
       cells       dict containing element vertex numbers
       bndDomain   dict id -> list[lower,upper] (or id -> str) where lower and upper describe the bounding box of the boundary section
       bndSegments dict id -> list of lists containing vertex numbers of boundary segments
       periodic    string containing periodic boundary transformation
       dim         dimension of grid
       computeFirstIndex if true, first vertex index is computed bases on all cells vertices (default is on)

    Returns
        String containing the mesh description in DGF format.
    """
    import numpy as np

    if dim is None:
        if "tetra" in cells or "hexa" in cells:
            dim = 3
        elif "quad" in cells or "triangle" in cells:
            dim = 2
        else:
            dim = len(points[0])

    simplex = "triangle" if "triangle" in cells else None
    if dim == 3 and "tetra" in cells:
        simplex = "tetra"

    # default numbering is 0 based
    firstVertexIndex = 0
    for key, cellVertices in cells.items():
        firstVertexIndex = min(firstVertexIndex, np.min(cellVertices))

    dgf = "DGF\nVertex\n"
    # if first vertex index is not 0 we need to add this here
    if firstVertexIndex > 0:
        dgf += f"firstindex {firstVertexIndex}\n"

    for p in points:
        for i in range(dim):
            dgf += str(p[i]) + " "
        dgf += "\n"
    dgf += "#\n\n"

    if simplex is not None:
        dgf += "Simplex\n"
        for t in cells[simplex]:
            for v in t:
                dgf += str(v) + " "
            dgf += "\n"
        dgf += "#\n\n"

    if "quad" in cells:
        dgf += "Cube\n"
        # gmsh has a different reference quadrilateral
        vxmap = [0, 1, 3, 2]  # flip vertex 2 and 3
        for t in cells["quad"]:
            for i in range(4):
                dgf += str(t[vxmap[i]]) + " "
            dgf += "\n"
        dgf += "#\n\n"

    # boundary segments
    if bndSegments is not None:
        assert isinstance(
            bndSegments, dict
        ), "Expecting a dictionary for boundary domain"
        dgf += "BoundarySegments\n"
        for bndid, bndsegs in bndSegments.items():
            for segment in bndsegs:
                dgf += str(bndid)
                for vx in segment:
                    dgf += " " + str(vx)
                dgf += "\n"
        dgf += "#\n\n"

    # boundary domain section
    if bndDomain is not None:
        assert isinstance(bndDomain, dict), "Expecting a dictionary for boundary domain"
        dgf += "BoundaryDomain\n"
        for bndid, bnd in bndDomain.items():
            if isinstance(bnd, str):
                if bnd == "default":
                    dgf += bnd + " " + str(bndid) + "\n"
                else:
                    dgf += str(bndid) + " " + bnd + "\n"
            else:  # tuple or list
                # has to be either list or tuple here
                assert isinstance(bnd, (tuple, list))
                assert len(bnd) == 2  # a lower left and upper right corner
                dgf += str(bndid)
                for coord in bnd:
                    assert len(coord) == dim  # should a coordinate in the domain
                    for c in coord:
                        dgf += " " + str(c)
                dgf += "\n"

        dgf += "#\n"

    # periodic boundaries
    if periodic is not None:
        dgf += periodic

    return dgf


# a simple projection for boundary ids
def projectBoundaryIds(gridView):
    """
    Parameter:
       gridView     a grid view to project boundary ids for

    Returns:
       a piecewise discrete function containing values corresponding to
       adjacent boundaries. `0` refers to interior elements.
    """
    import io, sys
    from dune.fem.space import finiteVolume
    from dune.generator import algorithm

    code = """
    #include <dune/fem/misc/boundaryidprovider.hh>

    template <class GridView, class Intersection>
    int boundaryId( const GridView& gv, const Intersection& i )
    {
      return Dune::Fem::boundaryId( gv, i );
    }
    """

    space = finiteVolume(gridView)
    bndFunc = space.interpolate([0], name="bndId")

    bndId = None

    for e in gridView.elements:
        idx = gridView.indexSet.index(e)
        if e.hasBoundaryIntersections:
            for i in gridView.intersections(e):
                if bndId is None:
                    bndId = algorithm.load("boundaryId", io.StringIO(code), gridView, i)
                if i.boundary:
                    # hack: assumption is that index and dof correspond
                    # which is only the case in this specific situation
                    bndFunc.dofVector[idx] = bndId(gridView, i)

    return bndFunc
