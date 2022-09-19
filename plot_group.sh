#!/bin/bash

#setting the parameters
stepLimit=1e6
stepDuration=1e-5
skipSteps=0
measurementInterval=1e-1
nGreenParticles=2000  
nRedParticles=2000
areaFractionList=(0.7 0.8 0.82 0.85 0.88)
redD=3
greenD=3
greenPersistentD=3
k=100
tau=0.05
PeList=(80 120 150)
potentialRange=1
LennardJones=1
turnAround=1
redRedAdhesionMult=0
greenGreenAdhesionMutl=0
redGreenAdhesionMult=0

TARGET_FOLDER="/home/marius/PhD/CellMotility/agent_simulation/new_output/turnAround_persistence"

#Loop over parameter values
for Pe in "${PeList[@]}"
do 
for areaFraction in "${areaFractionList[@]}"
do
filepath="${TARGET_FOLDER}/areaFraction_${areaFraction}_Pe_${Pe}"


# python3 /home/marius/PhD/CellMotility/analysis/animation_simulation.py "${filepath}"
python3 /home/marius/PhD/CellMotility/analysis/final_snapshot_simulation.py "${filepath}"
# python3 /home/marius/PhD/CellMotility/analysis/mixing_index_simulation.py "${filepath}"
# python3 /home/marius/PhD/CellMotility/analysis/cluster_analysis_simulation.py "${filepath}"

done
done



#Plot the mixing index
# python3 ../analysis/plot_mixing_index_simulation.py "${TARGET_FOLDER}/mixingIndex"


