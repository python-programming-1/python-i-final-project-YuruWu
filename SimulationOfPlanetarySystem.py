import numpy as np
import math
import matplotlib.pyplot as plt
def hamiltonian_solve(dq_dt, dp_dt, q_0, p_0, h, N):
    """ Solves for dynamics of Hamiltonian system

    Args:
        dq_dt: Partial derivative of coordinates 
        dp_dt: Derivative of momenta 

    Kwargs:
        q_0: Initial position (float for d=1, ndarray for d>1) set to 0.0 as default
        p_0: Initial momentum (float for d=1, ndarray for d>1) set to 1.0 as default
        h: Step size (float) set to 0.1 as default
        N: Number of steps (int) set to 100 as default

    Returns:
        Q: Numpy array of positions
        P: Numpy array of momenta
    """
    
    P=np.zeros((N+1,2))
    Q=np.zeros((N+1,2))

    Q[0] =q_0
    P[0] =p_0
    for n in range(N):
        P[n+1]=P[n]+h*dp_dt(Q[n],P[n])
        Q[n+1]=Q[n]+h*dq_dt(Q[n],P[n+1])
        
    
    return  Q


G = 6.67408e-11
M = 1.989e30
m = 5.972e24


def dp_dt(q,p):
  modulus=0
  for i in q:
    modulus=modulus+i*i
  modulus=math.sqrt(modulus)
  return -7.9277e44/math.pow(modulus,3)*q

def dq_dt(q,p):
 
    return p/5.97198e24

result=hamiltonian_solve(dq_dt,dp_dt,q_0=np.array([0.,1.471e11]),p_0=np.array([1.809e29,0.]),h=10000,N=3500000)

x_coordinate_earth=[]
y_coordinate_earth=[]
for i in result:
  x_coordinate_earth.append(i[0])
for i in result:
  y_coordinate_earth.append(i[1])

print(max(x_coordinate_earth),max(y_coordinate_earth))
plt.figure()
plt.plot(x_coordinate_earth,y_coordinate_earth)
plt.show()
