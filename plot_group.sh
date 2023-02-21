#!/bin/bash

#setting the parameters
stepLimit=1e6
stepDuration=1e-5
skipSteps=0
measurementInterval=1e0
nGreenParticles=2000  
nRedParticles=2000
areaFractionList=(0.5)
redD=3
greenD=3
greenPersistentD=0.1
k=50
tau=0.1
PeList=(40 80 120 160 200)
potentialRange=1 #1.10868
LennardJones=1
turnAround=1
redRedAdhesionMult=0
greenGreenAdhesionMutl=0
redGreenAdhesionMult=0

TARGET_FOLDER="/home/marius/PhD/CellMotility/agent_simulation/output_23_01/dCIL_persistence_initial_blob"

#Loop over parameter values
for Pe in "${PeList[@]}"
do 
for areaFraction in "${areaFractionList[@]}"
do
filepath="${TARGET_FOLDER}/A_${areaFraction}_Pe_${Pe}"

python3 ../analysis/write_mixing_index_sim.py $filepath
python3 ../analysis/write_clustering_sim.py $filepath
# python3 ../analysis/plot_last_frame_sim.py $filepath
# python3 ../analysis/animation_sim.py $parameterFile

done
done



#Plot the mixing index
# python3 ../analysis/plot_mixing_index_sim.py "${TARGET_FOLDER}/mixingIndex"


