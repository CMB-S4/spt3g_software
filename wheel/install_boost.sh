#!/bin/bash

set -e

# Location of this script
pushd $(dirname $0) >/dev/null 2>&1
scriptdir=$(pwd)
popd >/dev/null 2>&1
echo "Wheel script directory = ${scriptdir}"

PREFIX=${scriptdir}/deps

boost_version=1_86_0
boost_dir=boost_${boost_version}
boost_pkg=${boost_dir}.tar.bz2

if [ ! -e ${boost_pkg} ]; then
    echo "Fetching boost..."
    curl -SL "https://archives.boost.io/release/1.86.0/source/${boost_pkg}" -o "${boost_pkg}"
fi

if [ ! -e ${boost_dir} ]; then
    echo "Unpacking boost..."
    tar xjf ${boost_pkg}
fi

pyincl=$(for d in $(python3-config --includes | sed -e 's/-I//g'); do echo "include=${d}"; done | xargs)

echo "Building boost..."
cd ${boost_dir}
./bootstrap.sh --prefix=${PREFIX} --with-python=$(which python3) --with-python-root=$(python3-config --prefix)
./b2 -j2 --layout=tagged ${pyincl} variant=release threading=multi link=shared runtime-link=shared install