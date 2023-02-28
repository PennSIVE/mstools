#!/bin/bash
# run this after building
arch=$(uname -m | sed 's/x86_64/amd64/g') 
docker run --rm -v $PWD:/work pennsive/mstools:4.2-$arch bash /work/dev/run_pkg_tests.sh