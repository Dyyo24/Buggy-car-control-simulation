from BuggySimulator import *

# initial your buggy
def initail(traj):
	v = vehicle(vehicle.state(X=traj[0,0],
	                        Y=traj[0,1], 
	                        ))
	currentState = v.state
	return v, currentState