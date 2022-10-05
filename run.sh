#!/bin/bash

#Crash when any command fails
set -e

#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
if [ $# == 1 ]; then
    parameterFile="$1"
else 
    parameterFile="/home/marius/PhD/CellMotility/agent_simulation/new_output/adhesion/areaFraction_0.8_Pe_160"
fi

#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build --config Debug

#Execute
# nReps=3
# if [ $nReps == 0 ] 
# then
#     ./build/src/agent_simulation $parameterFile "_tracks.csv" #2>>"${parameterFile}"
# else
#     for i in $(seq $nReps)
#     do
#         ./build/src/agent_simulation $parameterFile "_tracks_${i}.csv" #2>>"${parameterFile}"
#     done
# fi

# Ananlyse tracks
python3 ../analysis/cluster_analysis_simulation.py $parameterFile
python3 ../analysis/mixing_index_simulation.py $parameterFile

#Analysis
# echo "Plotting the tracks..."
# python3 ../analysis/final_snapshot_simulation.py $parameterFile
# python3 ../analysis/animation_simulation.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_simulation.py $parameterFile