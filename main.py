from BuggySimulator import *
import numpy as np
from controller import *
from util import *
import matplotlib.pyplot as plt
from Evaluation import *


# get the trajectory
traj = get_trajectory('buggyTrace.csv')
# initial the Buggy
vehicle = initail(traj,0)
n = 15000
X = []
Y = []
delta = []
xd = []
yd = []
phi = []
phid = []
deltad = []
F = []
minDist =[]
e1 = []
e2 = []

k = []
e1 = []
e1d = []
e2 = []
e2d = []

'''
your code starts here
'''

# preprocess the trajectory
passMiddlePoint = False
nearGoal = False
for i in range(n):
    command, error = controller(vehicle)
    vehicle.update(command = command)
    if i % 80 == 0:
        print('X',vehicle.state.X)
        print('Y',vehicle.state.Y)
        _,minIdx = find_nearest_points(vehicle.state.X, vehicle.state.Y, traj)
        print('Traj',traj[minIdx])
        print(' ')

    # termination check
    disError,nearIdx = find_nearest_points(vehicle.state.X, vehicle.state.Y, traj)
    stepToMiddle = nearIdx - len(traj)/2.0
    if abs(stepToMiddle) < 100.0:
        passMiddlePoint = True
        print('middle point passed')
    nearGoal = nearIdx >= len(traj)-50
    if nearGoal and passMiddlePoint:
        print('destination reached!')
        break
    # record states

    e1.append(error[0])
    e1d.append(error[1])
    e2.append(error[2])
    e2d.append(error[3])
    
    X.append(vehicle.state.X)
    Y.append(vehicle.state.Y)
    delta.append(vehicle.state.delta)
    xd.append(vehicle.state.xd)
    yd.append(vehicle.state.yd)
    phid.append(vehicle.state.phid);
    phi.append(vehicle.state.phi)
    deltad.append(command.deltad)
    F.append(command.F)
    minDist.append(disError)
    k.append(vehicle.state.k/1000)
#output = np.array([xd,yd,phid,delta,X,Y,phi]).T
#np.save("24-677_Project_2_BuggyStates_Dyyo.npy",output)
    
showResult(traj,X,Y,delta,xd,yd,F,phi,phid,minDist)
showerror(traj,X,Y,e1,e1d,e2,e2d)
evaluation(minDist, traj, X, Y,3)