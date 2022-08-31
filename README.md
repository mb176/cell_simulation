# Quick Start

## Runing locally

- Create a parameter file: This file should have no file type and be of the form (alignment is optional, see notes below on the exact formatting):
<pre>
stepLimit               : 1e5 
nGreenParticles         : 100 
nRedParticles           : 100
stepDuration            : 1e-5
measurementInterval     : 1e-1
areaFraction            : 0.5
redD                    : 2
greenD                  : 2
greenPlusD              : 0.01
k                       : 100
tau                     : 2
Pe                      : 70
sigma                   : 1.0
LennardJones            : 1
skippedSteps            : 0
</pre>

- Run run.sh and either give the parameter file directly as command-line input, or set it as "parameterFile" inside the file. The simulation will output a file called "\<parameterFile>_tracks.csv". You can use the python scripts from analysis to analyse and visualise the data, this is done by simply handing over the parameter file to the python script, which will then read the parameters and find the trackfile by its name.

- For larger simulation use run_group.sh: You can select an array of values and it will create appropriatly named paramter files and run them. You can also tell it to run python scripts as described above.

## Running on the math cluster

- You have to ssh into the math cluster and set up a folder containing all the contents of agent_simulation and analysis (except the output of previous simulations).

- To run a single file modify run_cluster.sh to point to the correct parameter file and specify the job and then submit it by calling
> qsub run_cluster.sh

- For larger simulations use submit_batch.sh: Similarly to run_group.sh it allows you to generate paramter files for an array of values and then submit them all to the queue.

# Parameters

The following parameters are given in the parameter file to set up the simulation:

- stepLimit: The number of time steps the simulation performs. Together with stepDuration that defines the duration of the simulation.

- nGreenParticles: Number of green particles in the system. Green particles change persistence and direction when touching a red particle. For the period of heightened persistence they are decpicted as purple. 

- nRedParticles: Number of red particles. These particles don't undergo persistence change, but they do change direction upon contact with green particles. All particles are initialised on random positions. Overlaps are dealt with by minimising the repulsive potential energy via gradient descent.

- stepDuration: The duration of a single step in the simulation in the natural units of the simulation, r0^2/D. This means one time unit is roughly the time it takes for the particle to pass a particle diameter by diffusion, i.e. without any persistence.

-  measurementInterval: The time between measurements in natural units. Here a measurement refers to a timestep where the positions and states of all the particles are recorded and written to the parameter file. This has to be a multiple of the step duration, so that there are time steps that end exactly on the measurement times.

- areaFraction: The fraction of the simulation space that is covered by particles. The particles are sphererical, have diameter one in our units and for this value we assume no overlap.

- redD: The angular diffusion constant for red particles, expressed in natural units. Higher values mean lower persistence.

- greenD: The angular diffusion constant for unexcited green particles.

- greenPlusD: The angular diffusion constant for green particles once they reach an excited state by touching a red particle.

- k: The spring constant when a harmonic potential is used. If a Lennard-Jones potential is used, the force is multiplied by this value.

- tau: The duration of the excited state, which green particles enter upon contact with red particles. After this time they revert back to unexcited green particles.

- LennardJones: Boolean value that decides whether or not to use a Lennard-Jones interaction between particles. If this is zero the simulation uses a much softer harmonic potential.

- skippedSteps: This sets the number of steps for which no position measurements will be taken. This can be used to skip the initial period, where the system settles into a steady state, to avoid unecessaryly larege track files.

After the start of the simulation the program will add some further parameters to the list at run time:

- nParticles: The total number of particles in the system.
- Length: The side-length of the system in natural units. This is determined based on areaFraction and nParticles.
- Seed: The seed to used initialise the random number generator, for reproduceability.

Finally there are some hidden parameters, which are determined in agent_simulation_config.h. These are set to sensible values which usually don't need to be changed:

- SEED: Here you can choose the seed for initialising the random number generator.
- STEP_SIZE: This refers to the initial step size of the GSL minimizer which we use to reduce the overlap of the initial positions.
- LINE_TOL: The error tolerance for the overall GSL minimisation.
- TOL: Error tolerance for the overall GSL minisation.
- MAX_ITER: Maximum number of iterations we use for the inital minimisation.
- DEBUG: If defined the simulation is run in debug mode, which means more verbose output and extra assert cases.
- MAX_STEP_DISPLACEMENT: This is used in debug mode to check that the particle displacement is never larger than this value. This is a useful check to make sure that our step duration is small enough that particles don't land on top of each other and then bounce around uncontrolable due to the high force close to the center of the potential.


# Further notes

- The area on which the simulation is run always has the origin in the bottom left corner, i.e. extends into the top right quadrant

- All lengths are multiples of the cell diameter, i.e. the area covered by a cell is 0.785398.

- When running the code will append your parameter file, first adding additional parameters that are chosen at run run time and then adding a log of the computation. You can also have your error messages redirected there.

- About the parameter file: It cannot contain empty lines, and there has to be exactly one white space between : and the value. The maximum line length is 200 characters, and each variable can have at most 50 characters. The ordering matters and has to match the one used by SetParameters().
