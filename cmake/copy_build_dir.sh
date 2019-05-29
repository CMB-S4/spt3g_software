#!/bin/sh

export build_dir=${1}

rsync -mavL $build_dir/ $build_dir/tarball --exclude=".*" --exclude=CMakeFiles --exclude=tarball --exclude=Makefile --exclude="*.cmake" --exclude="CMakeCache.txt" --exclude="Testing" --exclude="*.g3" --exclude="*.g3.gz" --exclude="*.pkl"
