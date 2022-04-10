#!/bin/bash
#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/test/test_parameters"

# gcc -lgsl -lgslcblas -lm main.c -o main   #2>>"${parameterFile}"
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
cd ../analysis
# echo "Plotting the tracks..."
# python3 animation_simulation.py $parameterFile
# echo "Analysing the tracks..."
# python3 agent_simulation_analysis.py $parameterFile