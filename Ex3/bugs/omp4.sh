#!/bin/bash
ulimit -s unlimited
export OMP_STACKSIZE="20 m"
./omp4.exe
