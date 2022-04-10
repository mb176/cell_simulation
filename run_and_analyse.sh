#!/bin/bash

# parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HighDensity1/HighDensity1"
# parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HighDensityControl1/HighDensityControl1"
# parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/LowDensity1/LowDensity1"
# parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/LowDensityControl1/LowDensityControl1"
parameterFile="/home/marius/PhD/CellMotility/agent_simulation/validation/MIPS/ABP_MIPS_3"
# echo "Start compilation..."
# gcc -ffast-math -lgsl -lgslcblas -lm main.c -o main   #2>>"${parameterFile}"
# #-ffast-math for faster flops
# #-lgsl -lgslcblas to use GSL library
# #-lm to link math library
# ./main $parameterFile 

#Analysis
cd ../analysis
echo "Animating trajectory..."
python3 RDF_simulation.py $parameterFile
# echo "Analysing the tracks..."
# python3 MSD.py $parameterFile
echo "Done."