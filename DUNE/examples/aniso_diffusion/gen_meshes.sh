#!/bin/env bash

this_dir=$( dirname -- "${BASH_SOURCE[0]}" )
REPO_ROOT=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" )/../../..)" &> /dev/null && pwd )
orig_geo_path=$(realpath $REPO_ROOT/Firedrake/examples/scripts/aniso_diffusion/aniso_diffusion_DeluzetNarski.geo)

old_str="h = 0.1"
for h in 0.1 0.05 0.025 0.0125 0.00625 0.003125 0.0015625 0.00078125; do
    new_str="h = $h"
    tmp_geo="$this_dir/square_h${h}.geo"
    msh_path="$this_dir/square_h${h}.msh"
    echo Generating "$msh_path"
    sed < $orig_geo_path -e "s/$old_str/$new_str/" > "$tmp_geo"
    \rm -f "$msh_path"
    gmsh -2 "$tmp_geo" -o "$msh_path"
    \rm -f "$tmp_geo"
done
echo Done