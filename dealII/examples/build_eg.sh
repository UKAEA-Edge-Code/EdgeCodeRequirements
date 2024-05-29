#!/bin/env bash

this_dir=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" ))" &> /dev/null && pwd )
cmake_opts=""

#==============================================================================
build_eg() {
    # Build
    eg_dir="$this_dir/$1"
    if [ -d "$eg_dir" ]; then
        cmake_cmd="cmake $cmake_opts -B build ."
        cd "$eg_dir" &> /dev/null
        eval "$cmake_cmd"
        cd - &> /dev/null
        cd "$eg_dir/build" &> /dev/null
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
    cmake_opts="-DDEAL_II_DIR=$dealii_dir"
}
#==============================================================================

set_cmake_options

build_eg "$1"