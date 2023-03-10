#!/bin/bash

set -e #Stop script when any command fails

# Choose parameters
stepLimit=1e7
stepDuration=1e-5
skipSteps=0
measurementInterval=5e0
nGreenParticles=2000
nRedParticles=2000
areaFractionList=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8)
redD=3
greenD=3
greenPersistentD=0.1
kList=(50)
tau=0.02
PeList=(240 280 320)
potentialRange=1  #1.10868
LennardJones=1
turnAround=1
redRedAdhesionMult=0
greenGreenAdhesionMutl=0
redGreenAdhesionMult=0

# Choose simulation options
sed -i "/#define DIFFERENTIAL_CIL/c\#define DIFFERENTIAL_CIL" src/agent_simulation_config.h
sed -i "/#define CIL_DELAY/c\#define CIL_DELAY -1.0" src/agent_simulation_config.h
sed -i "/#define STICKY_CONTACTS/c\// #define STICKY_CONTACTS" src/agent_simulation_config.h
sed -i "/#define turnAroundVariation/c\// #define turnAroundVariation M_PI" src/agent_simulation_config.h
sed -i "/#define CIL_COOLDOWN_DURATION/c\// #define CIL_COOLDOWN_DURATION 0.02" src/agent_simulation_config.h
sed -i "/#define NON_DIFFERENTIAL_PERSISTENCE/c\// #define NON_DIFFERENTIAL_PERSISTENCE" src/agent_simulation_config.h

# Choose initial conditions
sed -i "/#define INITIAL_BLOB/c\// #define INITIAL_BLOB" src/agent_simulation_config.h
sed -i "/#define ONLY_GREEN_CIL/c\// #define ONLY_GREEN_CIL" src/agent_simulation_config.h
sed -i "/#define INITIAL_PHASE_SEGREGATION/c\// #define INITIAL_PHASE_SEGREGATION" src/agent_simulation_config.h

# Choose number of repetitions
nReps=1


TARGET_FOLDER="/home/ma/m/mpb19/CellMotility/agent_simulation/output_23_03/phasediagram/t_100"

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
k=$(bc <<< "scale=5; 0.05*$Pe")

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

#Plot the mixing index
# python3 ../analysis/plot_mixing_index_simulation.py "${TARGET_FOLDER}/mixingIndex"


