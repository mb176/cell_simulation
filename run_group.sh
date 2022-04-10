#!/bin/bash

#setting the parameters
stepLimit=100000
nGreenParticles=1000
nRedParticles=1000
stepDuration=0.0005 #200*60*0.00005 = 0.6 seconds
measurementInterval=1
declare -a areaFractionList
areaFractionList=(0.1 0.3 0.5)
redD=10
greenD=10
greenPlusD=1
# k = 100 #Set as 5xPe 
tau=0.25
declare -a PeList
PeList=(1 5 10 15 20 25)
sigma=1.1

TARGET_FOLDER="/home/marius/PhD/CellMotility/agent_simulation/output/mixing_index_2"

# Compile the simulation 
gcc -ffast-math -lgsl -lgslcblas -lm main.c -o main   #2>>"${parameterFile}"
#-ffast-math for faster flops
#-lgsl -lgslcblas to use GSL library
#-lm to link math library

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
# echo "k                       : $(($Pe*5))" >> "$filepath" 
# echo "tau                     : $tau" >> "$filepath" 
# echo "Pe                      : $Pe" >> "$filepath" 
# echo "sigma                   : $sigma" >> "$filepath" 

# #Start simulation
# ./main $filepath 

# #Calculate mixing index
# python3 ../analysis/calculate_mixing_index_simulation.py $filepath

# Get final frame
python3 ../analysis/track_animation.py $filepath

done
done

#Plot the mixing index
python3 ../analysis/plot_mixing_index_simulation.py "${TARGET_FOLDER}/mixingIndex"


