#!/bin/bash

#Crash when any command fails
set -e

#parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output/HigherDensity/HigherDensity"
if [ $# == 1 ]; then
    parameterFile="$1"
else 
    parameterFile="/home/marius/PhD/CellMotility/agent_simulation/output_new_CIL/test/test"
fi

#Genrate make file using cmake with gprof flags for compiler and linker
cmake -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg . -B build

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

# Get analysis
gprof ./build/src/agent_simulation gmon.out > ${parameterFile}_analysis.txt
