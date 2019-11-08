#!/bin/sh

vers=$(git describe --always --tags --dirty 2>/dev/null)
sed 's/\$Version\$/\$Version: '$vers'\$/g' <&0
