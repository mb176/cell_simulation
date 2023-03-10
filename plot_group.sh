#!/bin/bash

#setting the parameters
stepLimit=1e6
stepDuration=1e-5
skipSteps=0
measurementInterval=1e0
nGreenParticles=2000
nRedParticles=2000
areaFractionList=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
redD=3
greenD=3
greenPersistentD=0.1
kList=(1)
tau=0.02
PeList=(0 40 80 120 160 200 240 280 320)
potentialRange=1  #1.10868
LennardJones=1
turnAround=1
redRedAdhesionMult=0
greenGreenAdhesionMutl=0
redGreenAdhesionMult=0

TARGET_FOLDER="/home/marius/PhD/CellMotility/agent_simulation/output_23_03/phasediagram/t_100"

#Loop over parameter values
for Pe in "${PeList[@]}"
do 
for areaFraction in "${areaFractionList[@]}"
do
for k in "${kList[@]}"
do
filepath="${TARGET_FOLDER}/A_${areaFraction}_Pe_${Pe}"
#"${TARGET_FOLDER}/A_${areaFraction}_Pe_${Pe}"

# python3 ../analysis/write_mixing_index_sim.py $filepath
# python3 ../analysis/write_clustering_sim.py $filepath
python3 ../analysis/plot_last_frame_sim.py $filepath
# python3 ../analysis/animation_sim.py $parameterFile

done
done
done



#Plot the mixing index
# python3 ../analysis/plot_mixing_index_sim.py "${TARGET_FOLDER}/mixingIndex"


