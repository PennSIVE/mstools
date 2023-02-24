#!/bin/bash
# runs test fof pennsive packages
cd $(dirname $0)
APPENDLOG=0

function run_test_script () {
    test_script=$1

    echo "Performing test on $test_script" 

    # Run test and redirect output (append or overwrite)
    if [ ! ${APPENDLOG} -eq 1 ]; then
        Rscript $test_script 1> /dev/null 2> pkg_tests.err
        APPENDLOG=1
    else
        Rscript $test_script 1> /dev/null 2>> pkg_tests.err
    fi
    
    if [ ! $? -eq 0 ]; then
        echo "Performing test on $test_script (FAILED)"
        exit $?

    else
        echo "Performing test on $test_script (SUCCEEDED)"
        exit 0
    fi
}

run_test_script dev/test_rtapas.R 