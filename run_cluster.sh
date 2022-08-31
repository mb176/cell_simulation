#!/bin/bash
#PBS -N sim
#PBS -m a
#PBS -q standard

#Note:  -N sets the name for the job
#       -m sets the mails you receive (beginning, end, abort)
#       -q selects the queue the job is submitted to (standard, jumbo, etc)

#Set directory, since queue looses memory of where you submit your script
cd /home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/


parameterFile="/home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/output/test/test_parameters"



#Execute
# ./build/src/agent_simulation $parameterFile 2>>"${parameterFile}_error"

#Analysis
echo "Plotting the tracks..."
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/final_snapshot_simulation.py "${parameterFile}"
# python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/calculate_mixing_index_simulation.py "${parameterFile}"
python3 /home/clustor2/ma/m/mpb19/CellMotility/analysis/cluster_analysis_simulation.py "${parameterFile}"
