#!/bin/bash

#setting the parameters
stepLimit=1e6
nGreenParticles=2000
nRedParticles=2000
stepDuration=1e-5 #200*60*0.00005 = 0.6 seconds
measurementInterval=1e-1
declare -a areaFractionList
areaFractionList=(0.5)
redD=3
greenD=3
greenPlusD=0.3
k=100 #Set as 5xPe 
tau=0.25
declare -a PeList
PeList=(80 120)
sigma=1.0
LennardJones=1
skippedSteps=0

TARGET_FOLDER="/home/clustor2/ma/m/mpb19/CellMotility/agent_simulation/output/variedPersistence"

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
# echo $filepath
# echo "stepLimit               : $stepLimit" > "$filepath" 
# echo "nGreenParticles         : $nGreenParticles" >> "$filepath" 
# echo "nRedParticles           : $nRedParticles" >> "$filepath" 
# echo "stepDuration            : $stepDuration" >> "$filepath" 
# echo "measurementInterval     : $measurementInterval" >> "$filepath" 
# echo "areaFraction            : $areaFraction" >> "$filepath" 
# echo "redD                    : $redD" >> "$filepath" 
# echo "greenD                  : $greenD" >> "$filepath"
# echo "greenPlusD              : $greenPlusD" >> "$filepath" 
# echo "k                       : $k" >> "$filepath" #(($Pe*5))
# echo "tau                     : $tau" >> "$filepath" 
# echo "Pe                      : $Pe" >> "$filepath" 
# echo "sigma                   : $sigma" >> "$filepath" 
# echo "LennardJones            : $LennardJones" >> "$filepath" 
# echo "skippedSteps            : $skippedSteps" >> "$filepath"

#Modify the submission script (run_cluster.sh)
sed -i "/#PBS -N */c\#PBS -N ${filepath}" run_cluster.sh #-i for inplace, searches patter /.../ and changes it (c) to \...
sed -i "/parameterFile=/c\parameterFile=${filepath}" run_cluster.sh

#submit the job
qsub run_cluster.sh


done
done

#Plot the mixing index
# python3 ../analysis/plot_mixing_index_simulation.py "${TARGET_FOLDER}/mixingIndex"


