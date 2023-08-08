#!/bin/bash

#Crash when any command fails
set -e

#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
if [ $# == 1 ]; then
    parameterFile="$1"
else 
    parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output_23_07/comparison_clean_baseline/CIL_0.8"
fi

#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build

#Execute
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



# python3 /home/marius/PhD/CellMotility/analysis/animation_sim.py "${parameterFile}"
# vlc "${parameterFile}.mov"

python3 /home/marius/PhD/CellMotility/analysis/plot_last_frame_sim.py $parameterFile
python3 /home/marius/PhD/CellMotility/analysis/write_mixing_index_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/write_clustering_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_sim.py $parameterFile

# python3 /home/marius/PhD/CellMotility/analysis/plot_clustering_over_time_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_simulation.py $parameterFi


# # Copy this script into the folder
# cp /home/ma/m/mpb19/CellMotility/agent_simulation/submit_batch.sh $TARGET_FOLDER/submit_batch.sh

# #Handy commands
# # for i in {1830904..1830975}; do qdel $i; done
# # for i in {1830904..1830975}; do canceljob $i; done