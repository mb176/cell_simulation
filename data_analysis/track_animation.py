import matplotlib.pyplot as plt 
import numpy as np
import csv
import os
import matplotlib.animation as animation
from celluloid import Camera
from matplotlib.collections import PatchCollection
#python celluloid
#You can use np.arrays with named vectors, just like pandas!
# rain_drops = np.zeros(n_drops, dtype=[('position', float, (2,)),
#                                       ('size',     float),
#                                       ('growth',   float),
#                                       ('color',    float, (4,))])

# # Initialize the raindrops in random positions and with
# # random growth rates.
# rain_drops['position'] = np.random.uniform(0, 1, (n_drops, 2))

def readTracksFile(filename):
    measurementTimes = []
    positionMeasurements = []
    colorMeasurements = []
    with open(filename, mode="r") as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            measurementTimes.append(float(lines[0])) 
            positionMeasurements.append([])
            colorMeasurements.append([])
            nParticles = int((len(lines)-1)/3)
            for particleIdx in range(nParticles):
                colorMeasurements[-1].append(lines[1+3*particleIdx])
                positionMeasurements[-1].append([float(lines[1+3*particleIdx+1]), #x value
                                              float(lines[1+3*particleIdx+2])])   #y value

    return measurementTimes, colorMeasurements, np.array(positionMeasurements)
     
def plotColors(colorMeasurement): ##ToDo: Get right colors
    c = []
    for color in colorMeasurement:
        if color==" red":
            c.append('red')
        elif color==" green":
            c.append('green')
        elif color==" greenPlus":
            c.append('purple')
        else:
            assert False, "Invalid color value"
    return c

# 



SOURCE_FOLDER = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', "output"))+"/" #wtf python
NAME = "tracks"
measurementTimes, colorMeasurements, positionMeasurements = readTracksFile(SOURCE_FOLDER+NAME+".csv")
positionMeasurements = np.array(positionMeasurements)

xLength = 6 #ToDo: Get these from the parameter file
yLength = 6
R = 0.5
fig,ax = plt.subplots()
plt.axis([0, xLength, 0, yLength])
plt.axis("equal")

camera = Camera(fig)

for frameIdx in range(len(measurementTimes)):
    c = plotColors(colorMeasurements[frameIdx])
    for particleIdx in range(len(positionMeasurements[0])):
        circle = plt.Circle((positionMeasurements[frameIdx][particleIdx][0],
                                    positionMeasurements[frameIdx][particleIdx][1]), 
                                    radius=R, linewidth=0, color = c[particleIdx])
        plt.gca().add_artist(circle)
    camera.snap()
    #Delete artist before next frame
    for artist in plt.gca().lines + plt.gca().collections:
        artist.remove()
animation = camera.animate()
plt.show()
animation.save(SOURCE_FOLDER+"4particles.mov")


#AnimateTrajectories(measurementTimes, positionMeasurements, colorMeasurements)



# def AnimateTrajectories(measurementTimes, positionMeasurements, colorMeasurements):
#     xLength = 20
#     yLength = 10
#     R = 0.5
#     nFrames = len(measurementTimes)

#     fig, ax = plt.subplots()
#     ax.set_xlim(0, xLength)#, ax.set_xticks([]) 
#     ax.set_ylim(0, yLength)#, ax.set_yticks([])
#     ax.axis('equal')
    
#     x = positionMeasurements[0].transpose()[0]
#     y = positionMeasurements[0].transpose()[1]
#     c = plotColors(colorMeasurements[0])
#     scatter = ax.scatter(x, y, c=c) #ToDo: get actual cell radius
#     # circles = [plt.Circle((position[0],position[1]), radius=R, linewidth=0) for position in positionMeasurements[0]]
#     # scatter = PatchCollection(circles)
#     # ax.add_collection(scatter)
    

#     # scatter.set_offsets(positionMeasurements[1])
#     # plt.show()
#     ani = animation.FuncAnimation(fig, update_plot, frames=range(nFrames),
#                                 fargs=(scatter, positionMeasurements, colorMeasurements, measurementTimes,R))
    
#     #ani.save(SOURCE_FOLDER+"test"+".mov",fps=10)
#     plt.show()

# def update_plot(i, scatter, positionMeasurements, colorMeasurements, measurementTimes, R):
#         c = plotColors(colorMeasurements[i])
#         scatter.set_offsets(positionMeasurements[i])
#         scatter.set_array(c)
#         # circles = [plt.Circle((position[0],position[1]), radius=R, linewidth=0) for position in positionMeasurements[i]]
#         # scatter = PatchCollection(circles)        
#         #return scatter -> This function should not return anything!