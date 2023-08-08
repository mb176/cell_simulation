#!/bin/bash

#Crash when any command fails
set -e



#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
if [ $# == 1 ]; then
    parameterFile="$1"
else 
    parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output_23_07/realistic_parameters_CIL_distribution/movies/movie_script"
fi

# Choose simulation options
sed -i "/#define DIFFERENTIAL_CIL/c\//#define DIFFERENTIAL_CIL" src/agent_simulation_config.h
sed -i "/#define DIFFUSION_STRENGTH/c\#define DIFFUSION_STRENGTH 0.0" src/agent_simulation_config.h 
sed -i "/#define CIL_DELAY/c\#define CIL_DELAY 0.01" src/agent_simulation_config.h
sed -i "/#define STICKY_CONTACTS/c\// #define STICKY_CONTACTS" src/agent_simulation_config.h
sed -i "/#define TURN_AROUND_VARIATION/c\#define TURN_AROUND_VARIATION" src/agent_simulation_config.h
sed -i "/#define CIL_COOLDOWN_DURATION/c\//#define CIL_COOLDOWN_DURATION -1" src/agent_simulation_config.h
sed -i "/#define NON_DIFFERENTIAL_PERSISTENCE/c\//#define NON_DIFFERENTIAL_PERSISTENCE" src/agent_simulation_config.h

# Choose initial conditions
sed -i "/#define INITIAL_BLOB/c\// #define INITIAL_BLOB" src/agent_simulation_config.h
sed -i "/#define ONLY_GREEN_CIL/c\// #define ONLY_GREEN_CIL" src/agent_simulation_config.h
sed -i "/#define INITIAL_PHASE_SEGREGATION/c\// #define INITIAL_PHASE_SEGREGATION" src/agent_simulation_config.h

#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build --config Debug

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



python3 /home/marius/PhD/CellMotility/analysis/animation_sim.py "${parameterFile}"
vlc "${parameterFile}.mov"

# python3 /home/marius/PhD/CellMotility/analysis/plot_last_frame_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/write_mixing_index_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/write_clustering_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_sim.py $parameterFile

# python3 /home/marius/PhD/CellMotility/analysis/plot_clustering_over_time_sim.py $parameterFile
# python3 /home/marius/PhD/CellMotility/analysis/plot_mixing_index_simulation.py $parameterFi