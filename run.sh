#!/bin/bash

#Crash when any command fails
set -e

#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
if [ $# == 1 ]; then
    parameterFile="$1"
else 
    parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output_new_CIL/test/test"
fi

#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build --config Debug

# #Execute
nReps=0
if [ $nReps == 0 ] 
then
    ./build/src/agent_simulation $parameterFile "_tracks.csv" #2>>"${parameterFile}"
else
    for i in $(seq $nReps)
    do
        ./build/src/agent_simulation $parameterFile "_tracks_${i}.csv" #2>>"${parameterFile}"
    done
fi

# Ananlyse tracks
# python3 ../analysis/plot_clustering_over_time_sim.py $parameterFile
# python3 ../analysis/plot_mixing_index_sim.py $parameterFile

#Analysis
# echo "Plotting the tracks..."
# python3 ../analysis/final_snapshot_simulation.py $parameterFile
python3 ../analysis/animation_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_simulation.py $parameterFi