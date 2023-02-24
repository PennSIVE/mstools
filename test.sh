#!/bin/bash
docker run --rm -v $PWD:/work pennsive/mstools:4.2 bash /work/dev/run_pkg_tests.sh