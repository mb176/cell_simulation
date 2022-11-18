#!/bin/bash

set -e #Stop script when any command fails

#setting the parameters
stepLimit=1e6
stepDuration=1e-5
skipSteps=0
measurementInterval=1e0
nGreenParticles=2000
nRedParticles=2000
areaFractionList=(0.2 0.4 0.6 0.7)
redD=9
greenD=9
greenPersistentD=1
k=50
tau=0.02
PeList=(20 40 80 120 160)
potentialRange=0 #1.10868
LennardJones=1
turnAround=0
redRedAdhesionMult=0
greenGreenAdhesionMutl=0
redGreenAdhesionMult=0

nReps=10 #How often the simulation is to be repeated

TARGET_FOLDER="/home/ma/m/mpb19/CellMotility/agent_simulation/output_delayed_persistence/persistence_9_1/tau_0.02"

#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build --config Debug #2>"${parameterFile}_error"

#Loop over parameter values
for Pe in "${PeList[@]}"
do 
for areaFraction in "${areaFractionList[@]}"
do
filepath="${TARGET_FOLDER}/areaFraction_${areaFraction}_Pe_${Pe}"

# # Write parameter file
# redRedAdhesionMult=$(bc <<< "scale=5; 0.37176*$Pe/$k")
# greenGreenAdhesionMutl=$(bc <<< "scale=5; 0.37176*$Pe/$k")
# redGreenAdhesionMult=$(bc <<< "scale=5; 0.37176*$Pe/$k")

echo $filepath
echo "stepLimit               : $stepLimit" > "$filepath" 
echo "stepDuration            : $stepDuration" >> "$filepath" 
echo "skipSteps               : $skipSteps" >> "$filepath"
echo "measurementInterval     : $measurementInterval" >> "$filepath" 
echo "nGreenParticles         : $nGreenParticles" >> "$filepath" 
echo "nRedParticles           : $nRedParticles" >> "$filepath" 
echo "areaFraction            : $areaFraction" >> "$filepath" 
echo "redD                    : $redD" >> "$filepath" 
echo "greenD                  : $greenD" >> "$filepath"
echo "greenPersistentD        : $greenPersistentD" >> "$filepath" 
echo "k                       : $k" >> "$filepath"
echo "tau                     : $tau" >> "$filepath" 
echo "Pe                      : $Pe" >> "$filepath" 
echo "potentialRange          : $potentialRange" >> "$filepath" 
echo "LennardJones            : $LennardJones" >> "$filepath" 
echo "turnAround              : $turnAround" >> "$filepath"
echo "redRedAdhesionMult      : $redRedAdhesionMult" >> "$filepath"
echo "greenGreenAdhesionMutl  : $greenGreenAdhesionMutl" >> "$filepath"
echo "redGreenAdhesionMult    : $redGreenAdhesionMult" >> "$filepath"


#Modify the submission script (run_cluster.sh)
sed -i "/#PBS -N */c\#PBS -N ${filepath}" run_cluster.sh #-i for inplace, searches patter /.../ and changes it (c) to \...
sed -i "/parameterFile=/c\parameterFile=${filepath}" run_cluster.sh
sed -i "/nReps = /c\nReps=${nReps}" run_cluster.sh

#submit the job
qsub run_cluster.sh



done
done

#Plot the mixing index
# python3 ../analysis/plot_mixing_index_simulation.py "${TARGET_FOLDER}/mixingIndex"


