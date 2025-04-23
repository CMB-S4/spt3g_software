#!/bin/sh

export build_dir=${1}

for f in bin spt3g lib env-shell.sh; do
    rsync -mavL $build_dir/$f $build_dir/tarball
done
