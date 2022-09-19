#!/bin/bash
#PBS -N sim
#PBS -m a
#PBS -q standard

#Note:  -N sets the name for the job
#       -m sets the mails you receive (beginning, end, abort)
#       -q selects the queue the job is submitted to (standard, jumbo, etc)

#Set directory, since queue looses memory of where you submit your script
cd /home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/

parameterFile="/home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/new_output/adhesion/areaFraction_0.8_Pe_120"

#Execute
./build/src/agent_simulation $parameterFile 2>>"${parameterFile}"

#Analysis
echo "Plotting the tracks..."
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/animation_simulation.py "${parameterFile}"
python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/final_snapshot_simulation.py "${parameterFile}"
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/calculate_mixing_index_simulation.py "${parameterFile}"
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/cluster_analysis_simulation.py "${parameterFile}"
