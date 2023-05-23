#!/bin/bash

set -e #Stop script when any command fails

# Choose parameters
stepLimit=1e6
stepDuration=1e0
skipSteps=0
measurementInterval=1e5
nGreenParticles=2000
nRedParticles=2000
areaFractionList=(0.8)
redD=2
greenD=2
greenPersistentD=0
kList=(0.001)
tau=64
PeList=(0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5 )
potentialRange=1.20978631683  #1.10868
LennardJones=0
turnAround=1
redRedAdhesionMult=1
greenGreenAdhesionMutl=1
redGreenAdhesionMult=1

# Choose simulation options
sed -i "/#define DIFFERENTIAL_CIL/c\//#define DIFFERENTIAL_CIL" src/agent_simulation_config.h
sed -i "/#define CIL_DELAY/c\#define CIL_DELAY -1.0" src/agent_simulation_config.h
sed -i "/#define STICKY_CONTACTS/c\// #define STICKY_CONTACTS" src/agent_simulation_config.h
sed -i "/#define turnAroundVariation/c\// #define turnAroundVariation M_PI" src/agent_simulation_config.h
sed -i "/#define CIL_COOLDOWN_DURATION/c\#define CIL_COOLDOWN_DURATION 0.02" src/agent_simulation_config.h
sed -i "/#define NON_DIFFERENTIAL_PERSISTENCE/c\#define NON_DIFFERENTIAL_PERSISTENCE" src/agent_simulation_config.h

# Choose initial conditions
sed -i "/#define INITIAL_BLOB/c\// #define INITIAL_BLOB" src/agent_simulation_config.h
sed -i "/#define ONLY_GREEN_CIL/c\// #define ONLY_GREEN_CIL" src/agent_simulation_config.h
sed -i "/#define INITIAL_PHASE_SEGREGATION/c\// #define INITIAL_PHASE_SEGREGATION" src/agent_simulation_config.h

# Choose number of repetitions
nReps=1


TARGET_FOLDER="/home/ma/m/mpb19/CellMotility/agent_simulation/output_23_05/recreate_RnT/k_0.001"

#Genrate make file using cmake
cmake . -B build

#compile
cmake --build build --config Debug #2>"${parameterFile}_error"

#Loop over parameter values
for Pe in "${PeList[@]}"
do 
for areaFraction in "${areaFractionList[@]}"
do
for k in "${kList[@]}"
do
filepath="${TARGET_FOLDER}/A_${areaFraction}_Pe_${Pe}"

# # Write parameter file
# redRedAdhesionMult=$(bc <<< "scale=5; 0.37176*$Pe/$k")
# greenGreenAdhesionMutl=$(bc <<< "scale=5; 0.37176*$Pe/$k")
# redGreenAdhesionMult=$(bc <<< "scale=5; 0.37176*$Pe/$k")
# k=$(bc <<< "scale=5; 2*$Pe")

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
sed -i "/nReps=/c\nReps=${nReps}" run_cluster.sh

#submit the job
qsub run_cluster.sh


done
done
done

# Plot the mixing index
# python3 ../analysis/plot_mixing_index_simulation.py "${TARGET_FOLDER}/mixingIndex"


