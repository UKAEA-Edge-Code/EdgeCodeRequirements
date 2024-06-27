#!/bin/env bash

this_dir=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" ))" &> /dev/null && pwd )
cmake_opts=""
build_mode="Release"

#==============================================================================
build_eg() {
    # Build
    eg_dir="$this_dir/$1"
    if [ -d "$eg_dir" ]; then
        build_dir="builds/$build_mode"
        cmake_cmd="cmake $cmake_opts -B $build_dir ."
        cd "$eg_dir" &> /dev/null
        eval "$cmake_cmd"
        cd - &> /dev/null
        cd "$eg_dir/$build_dir" &> /dev/null
        make
        #make run
        cd - &> /dev/null
    else
        echo "No example dir at $eg_dir"
        exit 2
    fi
}

set_cmake_options()
{
    # Get dealii location
    dealii_dir="$(spack location -i dealii)"
    if [ $? -ne 0 ]; then
        echo "Couldn't find 'dealii' spack pkg"
        exit 1
    fi
    cmake_opts="-DDEAL_II_DIR=$dealii_dir -DCMAKE_BUILD_TYPE=$build_mode"
}
#==============================================================================

eg_name="$1"

if [ $# -eq 2 ]; then
    build_mode="$2"
elif [ $# -eq 1 ]; then
    build_mode="Release"
else
    echo "usage: build_eg.sh [eg_dir] <build_mode>"
    exit 2
fi
set_cmake_options
build_eg "$eg_name" "$build_mode"