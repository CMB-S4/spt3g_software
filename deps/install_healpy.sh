#!/bin/bash
#
# Install a relocatable healpy into the current test environment.
#
# When PyPI has no binary healpy wheel for this Python (e.g. new releases),
# pip builds healpy from the sdist.  That build compiles its bundled
# cfitsio / libsharp / healpix_cxx as shared libraries in a temporary
# directory and links against them there, so the installed package can't
# find them at import time.  Here we build in a persistent directory and
# run auditwheel/delocate on the result, just like the official wheels.

set -e

PYTHON=${PYTHON:-python}

workdir=$(mktemp -d)
echo "healpy work directory = ${workdir}"

# Nothing to do if a binary wheel exists for this interpreter/platform.
if ${PYTHON} -m pip download --only-binary=:all: --no-deps healpy \
    -d "${workdir}/probe" >/dev/null 2>&1; then
    echo "Binary healpy wheel available; skipping source build."
    exit 0
fi

echo "No binary healpy wheel found; building from source..."

${PYTHON} -m pip download --no-binary=:all: --no-deps healpy -d "${workdir}"
tar xzf "${workdir}"/healpy-*.tar.gz -C "${workdir}"
srcdir=$(ls -d "${workdir}"/healpy-*/)

# Build from the unpacked tree so setuptools' build/ directory (which holds
# the bundled shared libraries) survives until the wheel is repaired.
${PYTHON} -m pip wheel --no-deps -w "${workdir}/dist" "${srcdir}"

# Directories containing the bundled shared libraries
libdirs=$(find "${srcdir}/build" \( -name 'lib*.so*' -o -name 'lib*.dylib' \) \
    -exec dirname {} \; | sort -u | tr '\n' ':')
echo "Bundled library directories = ${libdirs}"

mkdir -p "${workdir}/fixed"
case "$(uname -s)" in
    Linux)
        LD_LIBRARY_PATH="${libdirs}${LD_LIBRARY_PATH}" \
            auditwheel repair -w "${workdir}/fixed" "${workdir}"/dist/healpy-*.whl
        ;;
    Darwin)
        ${PYTHON} -m pip install delocate
        DYLD_LIBRARY_PATH="${libdirs}${DYLD_LIBRARY_PATH}" \
            delocate-wheel -v -w "${workdir}/fixed" "${workdir}"/dist/healpy-*.whl
        ;;
    *)
        cp "${workdir}"/dist/healpy-*.whl "${workdir}/fixed/"
        ;;
esac

${PYTHON} -m pip install "${workdir}"/fixed/healpy-*.whl
