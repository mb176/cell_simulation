#!/bin/bash
#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
if [ $# == 1 ]; then
    parameterFile="$1"
else 
    parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/test/test_parameters"
fi


# gcc -lgsl -lgslcblas -lm main.c -o main #2>>"${parameterFile}"
#-ffast-math for faster flops
#-lgsl -lgslcblas to use GSL library
#-lm to link math library


#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build --config Debug

#Execute
./build/src/agent_simulation $parameterFile #2>>"${parameterFile}"

#Analysis
echo "Plotting the tracks..."
# python3 ../analysis/final_snapshot_simulation.py $parameterFile
# python3 ../analysis/animation_simulation.py $parameterFile