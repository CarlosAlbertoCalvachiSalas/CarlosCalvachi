import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib as mpl
from tqdm import tqdm


Min, Max, N = 0.,40.,11
x = np.linspace(Min,Max,N)
y = x.copy()
h = x[1]-x[0]

def h1(y):
    return 100.
def h2(y):
    return 0.
def h3(x):
    return 0.
def h4(x):
    return 0.

def InitT():
    
    T = np.zeros((N,N))
    
    T[0,:] = h1(y)
    T[-1,:] = h2(y)
    
    T[:,0] = h3(x)
    T[:,-1] = h4(x)
    
    return T


def GetRelaxation(T, Nit = int(1e5), omega = 1.9, tolerancia = 1e-2):
    
    itmax = 0
    
    for it in range(Nit):
        
        dmax = 0.
        
        for i in range(1, len(x)-1):
            for j in range(1, len(y)-1):
                tmp = 0.25*(T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1])
                r = omega*(tmp - T[i,j])
                
                T[i,j] += r
                
                if np.abs(r) > dmax:
                    dmax = r
       # print(T)
       # print(it)
        
        if np.abs(dmax) < tolerancia:
            #print(it)
            itmax = it
            break
            
    return T,itmax

omegas = np.linspace(1, 1.9, 1000)

def getIterations(omega):
	T = InitT()
	Tf1, itmax =  GetRelaxation(T, omega = omega)

	return itmax

iterations = np.vectorize(getIterations)(omegas)

print('omega más óptimo:', np.average(omegas[iterations == np.min(iterations)]))

plt.plot(omegas, iterations, color = 'blue')
plt.title('Optimization of the Over-relaxation Method')
plt.grid()
plt.ylabel('Iterations')
plt.xlabel(r'$\omega$')
plt.show()


