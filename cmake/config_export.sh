#!/bin/sh

git config filter.exportversion.smudge cmake/smudge_version.sh
git config filter.exportversion.clean cmake/clean_version.sh
