'''
This is a realization of the controller from Vehicle Dynamics and Control by Rajesh Rajamani.
'''
from BuggySimulator import *
import numpy as np
import scipy
import control
from scipy.ndimage import gaussian_filter1d
from util import *

def dlqr(A,B,Q,R):

    X = np.matrix(scipy.linalg.solve_discrete_are(A, B, Q, R))
 
        #compute the LQR gain
    K = np.matrix(scipy.linalg.inv(B.T*X*B+R)*(B.T*X*A))
 
    eigVals, eigVecs = scipy.linalg.eig(A-B*K)
 
    return K, X, eigVals


def PID(p,i,d,target,current,dt):
    u = p*(target - current) + i*(target - current)*dt + d*(target - current)/dt
    return u



traj = get_trajectory('buggyTrace.csv')
lastF = 0

def curvature(point1,point2,point3,dt):
    
    cxprime = (point2[0] - point1[0])/dt
    cnext_xprime = (point3[0] - point2[0])/dt
    cxdprime = (cnext_xprime - cxprime)/dt
    
    cyprime = (point2[1] - point1[1])/dt
    cnext_yprime = (point3[1] - point2[1])/dt
    cydprime = (cnext_yprime - cyprime)/dt
    
    k = (cxprime*cydprime-cyprime*cxdprime)/(cxprime**2+cyprime**2)**1.5
    return k  
    

def controller(vehicle):
    
    ca = vehicle.Ca
    m = vehicle.m
    lr = vehicle.lr
    lf = vehicle.lf
    Iz = vehicle.Iz
    f = vehicle.f
    g = vehicle.g
    vx = vehicle.state.xd
    x = vehicle.state.X
    y = vehicle.state.Y
    yd = vehicle.state.yd
    phi = vehicle.state.phi
    phid = vehicle.state.phid
    delta = vehicle.state.delta

    'initialize A,B'
    A = np.array([[0,1,0,0],[0,-(4*ca)/(m*vx),4*ca/m,2*ca*(lr-lf)/m/vx],
                  [0,0,0,1],[0,-2*ca*(lf-lr)/Iz/vx,2*ca*(lf-lr)/Iz,-2*ca*(lf**2+lr**2)/Iz/vx]])
    B = np.array([[0],[2*ca/m],[0],[2*ca*lf/Iz]])
    
    'Get LQR'
  #  L = np.array([[2,0,0,0],[0,2,0,0],[0,0,2,0],[0,0,0,2]])
    R = 450
    Q = np.array([[30,0,0,0],[0,8,0,0],[0,0,2000,0],[0,0,0,8]])
    
    #Ad,Bd,Cd,Dd,_ = scipy.signal.cont2discrete((A,B,np.diag([1,1,1,1]),np.zeros(B.shape)),dt=0.05)
   # K,X,ev = dlqr(Ad,Bd,Q,R)
 
    P = scipy.linalg.solve_continuous_are(A, B, Q, R)
    K = -1/R*B.T@P
    'Calculate K the fastest way'
    # LQR return a K in some steps
    
    # calculate error

    
    distance, index = find_nearest_points(x,y,traj)
   # print("The tracking point is:,", index)
    #print('index',index)
    lookahead = 50
    front_point = np.roll(traj,-lookahead)[index]
    front_point1 = np.roll(traj,-lookahead-2)[index]
    front_point2 = np.roll(traj,-lookahead-4)[index]
    
    'Calculate phi_desired'
    phi_desired = np.arctan2([front_point1[1]-front_point[1]],[front_point1[0]-front_point[0]])[0]
    phi_desired = wrap2pi(phi_desired)
    #phi_desired = np.arctan2([front_point[1]],[front_point[0]])[0]
    #print(phi_desired)
    
    'Calculate curvature'
    trajc = scipy.ndimage.filters.gaussian_filter1d(traj,sigma=1)
    dt = 1
    clookahead = 1     # for calculating error
    cfront_point = np.roll(trajc,-clookahead)[index]
    cfront_point1 = np.roll(trajc,-clookahead-1)[index]
    cfront_point2 = np.roll(trajc,-clookahead-2)[index]
    
#    cxprime = (cfront_point1[0] - cfront_point[0])/dt
#    cnext_xprime = (cfront_point2[0] - cfront_point1[0])/dt
#    cxdprime = (cnext_xprime - cxprime)/dt
#    
#    cyprime = (cfront_point1[1] - cfront_point[1])/dt
#    cnext_yprime = (cfront_point2[1] - cfront_point1[1])/dt
#    cydprime = (cnext_yprime - cyprime)/dt
#    
#    k = (cxprime*cydprime-cyprime*cxdprime)/(cxprime**2+cyprime**2)**1.5
    k = curvature(cfront_point,cfront_point1,cfront_point2,dt)
    
    vehicle.state.k = k
    
    phid_desired = vx*k
    
    'Calculate e1'
    vecty = [front_point[0] - x, front_point[1] - y]    
#    vect_tan = [front_point1[0] - front_point[0], front_point1[1] - front_point[1]] 
    vect_perpen = [front_point1[1] - front_point[1], -front_point1[0] + front_point[0]]
#    print('vect_perpen',vect_perpen)
#    print('Vecty',vecty)
    # np.inner(vecty,vect_perpen.reshape(1,2))
#    e1 = ((front_point[0] - x)*(front_point1[1] - front_point[1]) + (front_point[1] - y)*(-front_point1[0] + front_point[0]))
#    e1 = e1/np.sqrt(vect_perpen[0]**2 + vect_perpen[1]**2)
    e1 = (vehicle.state.Y - traj[index][1] )/ np.cos(phi_desired)
    #e1 = (front_point[0] - x)*(front_point1[1] - front_point[1]) + (front_point[1] - y)*(-front_point1[0] + front_point[0])
    e2 = wrap2pi (phi - phi_desired)
    #print(e2)
    e1d = yd + vx *(e2)
    #e2d = wrap2pi(phid - phid_desired)
    e2d = phid - phid_desired
    #print(e2d)
    #print(vehicle.state.phi)
    state = np.array([e1,e1d,e2,e2d])
   # print("The error is ", state)
# =============================================================================
#     Fx = 1
#     if vehicle.state.x - traj[i][0]:
#         Fx = -Fx
# =============================================================================
    if state.shape == (4,1,1):
        state = state.squeeze().reshape((4,1))
    u = K@state
    
    #eigVals, eigVecs = scipy.linalg.eig(A+B*K)
    #print('eigVals', eigVals)
    
    dtime = 0.05
    
 #   deltadot = PID(0.4,0.3,-0.1,u[0],delta,dtime)
    deltadot = -(delta - u[0])/dtime            # how to calculate dt?
    #print('delta dot',wrap2pi(delta - u[0]))
    #print(deltadot)i
    'Return desired states'
    
#    
#    curve_lookahead = 100     # for calculating error
#    cfront_point1 = np.roll(traj,-clookahead)[index]
#    cfront_point2 = np.roll(traj,-clookahead-2)[index]
#    cfront_point3 = np.roll(traj,-clookahead-4)[index]
#    k1 = curvature(cfront_point1,cfront_point2,cfront_point3,dt)
#    
#    curve_lookahead_ref = 75
#    cfront_point1 = np.roll(traj,-curve_lookahead_ref)[index]
#    cfront_point2 = np.roll(traj,-curve_lookahead_ref-2)[index]
#    cfront_point3 = np.roll(traj,-curve_lookahead_ref-4)[index]
#    k2 = curvature(cfront_point1,cfront_point2,cfront_point3,dt)
#    
#    curve_lookahead_ref = 50
#    cfront_point1 = np.roll(traj,-curve_lookahead_ref)[index]
#    cfront_point2 = np.roll(traj,-curve_lookahead_ref-2)[index]
#    cfront_point3 = np.roll(traj,-curve_lookahead_ref-4)[index]
#    k3 = curvature(cfront_point1,cfront_point2,cfront_point3,dt)
    
    v = 20
#    if abs(abs(k2) - abs(k3)) > 3*abs(abs(k2) - abs(k1)):
#        v = 2.5
#        change = False
    
    #if e1 > 2.5:
    #    v = 2
#    if x > 240 and x < 370 and y < -165:
#        v = 25
    
    # limit speed on corners
    if x < 220 and x >0 and y > -125:         # first dash
        v = 30
    
    if x > 370 and x < 450 and y > -165 and y < -92:        # first curve
        v = 3.3
        
    if x >= 130 and x < 155 and y < -180:        # second curve
        v = 5
    if x > 30 and x < 130 and y < -180:         # turning at second curve
        v = 3.7

    if x <= 30 and x >= -35 and y < -100:
        v = 10
    if x < -35 and y < -60 and y > -150:    # third curve 
        v = 3.5
    if x >-80 and x < 0 and y > -100:
        v = 35                        # final dash 
    

#    F = lastF
#    if v < 0.1:
#        lastF = 0
        
    if vehicle.state.xd <= v:
        F = 6000
#       Ftarget = 6000
#       F = PID(0.3,0.3,-0.01,Ftarget,lastF,dt)
    else:
        F = -8000
#       Ftarget = -8000
#       F = PID(0.3,0.3,-0.01,Ftarget,lastF,dt)
#    lastF = F
#    if abs(phi) > 1.5:
 #       F = 0
 #       phi = 1.5
        
    command = vehicle.command(F, deltadot)
    
    return command,[e1,e1d,e2,e2d]
