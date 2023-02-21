#!/bin/bash
#PBS -N /home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/new_output/adhesion/areaFraction_0.3_Pe_40
#PBS -m a
#PBS -q standard

#Note:  -N sets the name for the job
#       -m sets the mails you receive (beginning, end, abort)
#       -q selects the queue the job is submitted to (standard, jumbo, etc)

#Crash when any command fails
set -e

#Set directory, since queue looses memory of where you submit your script
cd /home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/

parameterFile=/home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/new_output/adhesion/areaFraction_0.3_Pe_40


# Execute
nReps=1
if [ $nReps == 0 ] 
then
    ./build/src/agent_simulation $parameterFile "_tracks.csv" 2>>"${parameterFile}"
else
    for i in $(seq $nReps)
    do
        ./build/src/agent_simulation $parameterFile "_tracks_${i}.csv" 2>>"${parameterFile}"
    done
fi

# Ananlyse tracks
python3 ../analysis/write_mixing_index_sim.py $parameterFile
python3 ../analysis/write_clustering_sim.py $parameterFile


# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/animation_sim.py "${parameterFile}"
python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/plot_last_frame_sim.py "${parameterFile}"
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/calculate_mixing_index_simulation.py "${parameterFile}"
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/cluster_analysis_simulation.py "${parameterFile}"
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/write_histogram.py "${parameterFile}"
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/write_histogram_delta.py "${parameterFile}"