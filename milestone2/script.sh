#!/bin/bash

N=10000N
type=".txt"
for var in "8_threads_static"
do
    set -f
    a=($(echo $var | tr '_' '\n'))
    echo "${a[0]}"

    # echo "--- compiling serial code."
    # make MyAdvectionSerial
    # #filename="${var}${type}"
    filename="metrics.txt"
    # echo "--- running serial code for lax method."
    # touch ./logs/${N}/serial/lax/${filename}
    # ./my_advection_program_serial 1131 400 1.0 1.0e3 5.0e-7 2.85e-7 lax ${a[0]} > ./logs/${N}/serial/lax/${filename}
    # echo "--- running serial code for first_order method."
    # touch ./logs/${N}/serial/first_order/${filename}
    # ./my_advection_program_serial 1131 400 1.0 1.0e3 5.0e-7 2.85e-7 first_order ${a[0]} > ./logs/${N}/serial/first_order/${filename}
    # echo "--- running serial code for second_order method."
    # touch ./logs/${N}/serial/second_order/${filename}
    # ./my_advection_program_serial 1131 400 1.0 1.0e3 5.0e-7 2.85e-7 second_order ${a[0]} > ./logs/${N}/serial/second_order/${filename}
    # echo "\n"
    # echo "-------------------------------------------------------------------------------"
    # echo "\n"
    echo "--- compiling parallel code."
    make MyAdvection
    echo "--- running parallel code for lax method."
    #touch ./logs/${N}/parallel/lax/${filename}
    ./my_advection_program_parallel 10000 20000 1.0 1.0e6 5.0e-7 2.85e-7 lax ${a[0]} > ./logs/${N}/parallel/lax/${filename}
    echo "--- running parallel code for first_order method."
    #touch ./logs/${N}/parallel/first_order/${filename}
    # ./my_advection_program_parallel 10000 20000 1.0 1.0e6 5.0e-7 2.85e-7 first_order ${a[0]} > ./logs/${N}/parallel/first_order/${filename}
    # echo "--- running parallel code for second_order method."
    # #touch ./logs/${N}/parallel/second_order/${filename}
    # ./my_advection_program_parallel 10000 20000 1.0 1.0e6 5.0e-7 2.85e-7 second_order ${a[0]} > ./logs/${N}/parallel/second_order/${filename}
done