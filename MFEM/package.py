# Copyright 2013-2023 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install py-pymfem
#
# You can edit this file again by typing:
#
#     spack edit py-pymfem
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack.package import *


class PyPymfem(PythonPackage):
    """Python wrapper for MFEM"""

    homepage = "http://mfem.org/"
    pypi = "mfem/mfem-4.6.1.0.tar.gz"

    # FIXME: Add a list of GitHub accounts to
    # notify when the package is updated.
    # maintainers("github_user1", "github_user2")

    version("4.6.1.0", sha256="cf10d8162349dca1ac2f49e0f401a8786ab950a81bbff567ea542252763f11d2")

    depends_on("python@3.7:", type=("build", "run"))

    depends_on("py-setuptools", type="build")
    depends_on("swig@4.1.1:", type="build")
    depends_on("cmake", type="build")

    #depends_on("mfem")
    depends_on("py-numpy@1.20.0:", type=("build", "run"))
    depends_on("py-scipy", type=("build", "run"))
    depends_on("py-numba", type=("build", "run"))
    depends_on("py-mpi4py", type=("build", "run"), when="+mpi")
    #depends_on("py-numba-scipy", type=("build", "run"))

    # Can I get it to use preinstalled METIS, HYPRE, CUDA, gslib, suitesparse

    # Might require significant changes to be able to build against preinstalled MFEM, to reflect the fact the source directory isn't present. Only uses headers, but not sure if all of them get installed. More importantly, there is the problem of what versions are compatible with each other. It looks like it may be specific commits only, which would basically make it impossible.

    # Should probably patch so spack controls number of cores used to compile

    variant("mpi", default=True, sticky=True, description="Enable MPI parallelism")

    def global_options(self, spec, prefix):
        # FIXME: Add options to pass to setup.py
        # FIXME: If not needed, delete this function
        options = []
        return options

    def install_options(self, spec, prefix):
        if "+mpi" in spec:
            return ["--with-parallel"]
        return []
