#!/bin/bash

export build_dir=${1}

rsync -avL $build_dir/ $build_dir/tarball --exclude=".*" --exclude=CMakeFiles --exclude=tarball --exclude="*.cmake" --exclude="CMakeCache.txt"
