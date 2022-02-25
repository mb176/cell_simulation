#!/bin/bash
#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/test/test_parameters"
gcc -lm main.c -o main  #-lm to link math library
./main $parameterFile 