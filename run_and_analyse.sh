#!/bin/bash
#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/test/test_parameters"
gcc -lgsl -lgslcblas -lm main.c -o main  
#-lm to link math library
#-lgsl -lgslcblas to use GSL library
./main $parameterFile 

#Analysis
cd ../analysis
echo "Animating trajectory..."
python3 track_animation.py $parameterFile
echo "Analysing the tracks..."
python3 agent_simulation_analysis.py $parameterFile
echo "Done."