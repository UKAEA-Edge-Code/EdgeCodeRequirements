#!/bin/env bash

this_dir=$( cd -- "$(realpath $( dirname -- "${BASH_SOURCE[0]}" ))" &> /dev/null && pwd )
repo_url="git@github.com:dealii/dealii.git"
# Use 9.4 to match spack package
version="9.4"
tmp_dir="/tmp"
repo_name="dealii"
repo_dir="$tmp_dir/$repo_name"

cmake_opts=""
#==============================================================================
clone_repo(){
    # Already cloned?
    repo_exists=false
    if [ -d "$repo_dir" ]; then
        cd "$repo_dir" > /dev/null 2>&1
        if git tag > /dev/null 2>&1; then
            repo_exists=true
        fi
        cd - > /dev/null 2>&1
    fi

    # Clone if not
    if [ "$repo_exists" = false ]; then
        cd "$tmp_dir" > /dev/null 2>&1
        git clone -b dealii-$version "$repo_url" "$repo_name"
        cd - > /dev/null 2>&1
    fi
}

copy_and_build_eg() {
    # Copy
    eg_dir="$this_dir/step-$1"
    echo $eg_dir
    src_eg_dir="$repo_dir/examples/step-$1"
    if [ -d "$src_eg_dir" ]; then
        rsync_cmd="rsync -av --exclude='doc' $src_eg_dir $this_dir  &> /dev/null"
        eval "$rsync_cmd"
    fi
    # Build
    cmake_cmd="cmake $cmake_opts -B build ."
    cd "$eg_dir" &> /dev/null
    eval "$cmake_cmd"
    cd - &> /dev/null
    cd "$eg_dir/build" &> /dev/null
    make
    #make run
    cd - &> /dev/null
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

clone_repo

# Copy egs from the repo and build
for eg_num in 1 2 3 4 5; do
    copy_and_build_eg "$eg_num"
done