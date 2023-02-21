#!/bin/bash

#Crash when any command fails
set -e

#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
if [ $# == 1 ]; then
    parameterFile="$1"
else 
    parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output_23_02/Rosalbas_videos/test_no_bug"
fi

#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build --config Debug

# #Execute
nReps=1
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
python3 /home/marius/PhD/CellMotility/analysis/plot_last_frame_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_clustering_over_time_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_sim.py $parameterFile

python3 /home/marius/PhD/CellMotility/analysis/animation_sim.py "${parameterFile}"

# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_simulation.py $parameterFi
# python3 /home/marius/PhD/CellMotility/analysis/write_mixing_index_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/write_clustering_sim.py $parameterFile