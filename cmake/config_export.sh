#!/bin/sh

git config filter.exportversion.smudge "git describe --always --tags --dirty"
