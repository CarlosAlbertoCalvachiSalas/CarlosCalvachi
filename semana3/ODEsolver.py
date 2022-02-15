import numpy as np
import matplotlib.pyplot as plt
import sys

sys.setrecursionlimit(100000)

# Euler Rule

# Equations of the form y'(t) = f(t) given initial condition 
# y(0) = a, with a as a constant.


def ODE(y0, t0 = 0, tmax = 1, f=lambda t: 1, N = 10000):
	t = np.linspace(t0, tmax, N) # it includes t = t0 
	f = np.vectorize(f)	

	h = tmax/(N-1)

	ys = np.append([y0] , h*f(t[:N-1]))

	# There are two approaches. The recursive and array like
	# Array like

	#preArray = np.zeros(2*N - 1) 
	#repeats  = np.zeros(2*N - 1) 

	# N - 1 white spaces and N steps

	#preArray[::2]  	= ys
	#repeats[::2] 	= (np.arange(N) + 1)[::-1]
	#repeats[1::2]	= (np.arange(N-1) + 1)

	#array = np.repeat(preArray, np.int_(repeats)).reshape(N,N)

	#y = np.sum(array, axis = 0)

	y = np.cumsum(ys)

	plt.plot(t,y, color = 'blue')
	


#ODE(-1, t0 = -5, f=lambda t: np.sin(2*t**2))
#plt.show()

# Equations of the form f'(t) = f(t,y) can be solved this way


"""
def eulerODE(f, y0, t0 = 0, tmax = 1, N = 10000):
	t = np.linspace(t0, tmax, N) # it includes t = t0 
	h = tmax/(N-1)
	y = []

	def eulerRecursion(n):
		if(n == 0):
			y.append(y0)
			return y0
		else:
			yn_1 = eulerRecursion(n-1)
			yn = yn_1 + h*f(t[n-1],yn_1)
			y.append(yn)
			return yn
		
	eulerRecursion(N-1)

	y = np.array(y)

	return t, y


eulerODE(lambda t,y: y, 1) 

plt.show()
"""

def rungeKuttaODE(f, y0, t0 = 0, tmax = 1, N = 10000):
	t = np.linspace(t0, tmax, N) # it includes t = t0 
	h = tmax/(N-1)
	y = []

	def rungeKuttaRecursion(n):
		if(n == 0):
			y.append(y0)
			return y0
		else:

			yn_1 = rungeKuttaRecursion(n-1)
			k1 = h*f(t[n-1], yn_1)
			k2 = h*f(t[n-1]+(h/2), yn_1 + (k1/2))
			k3 = h*f(t[n-1]+(h/2), yn_1 + (k2/2))
			k4 = h*f(t[n-1]+h, yn_1 + k3)

			yn = yn_1 + ((k1+2*k2+2*k3+k4)/6)

			y.append(yn)
			return yn

	rungeKuttaRecursion(N-1)

	y = np.array(y)

	return t, y


#rungeKuttaODE(lambda t,y: y, 1) 

#plt.show()


k = 389.6
A = 0.01 
l = 0.3
R = 8.31446261815324 # NIST value
cv = (3/2)*R

C = (k*A)/(cv*l)


def heatTransfer(C = C):
	M = np.array([[-C, C], [C, -C]])
	return lambda t, y: M@y

t, yArray = rungeKuttaODE(f = heatTransfer(), y0 = np.array([400,200]), t0 = 0, tmax = 4, N = 10000)


plt.plot(t, yArray[:, 0], label = r'$T_1$', color = 'blue')
plt.plot(t, yArray[:, 1], label = r'$T_2$', color = 'green')
plt.xlabel('t(s)')
plt.ylabel('T(K)')
plt.legend()
plt.show()


# Runge - Kutta algorithm

"""
For solving systems of equations what is needed is to
extrapolate f to give a vector to a vector. Then y0 is a vector

Lets take a harmonic oscillator

ø'' + omega*ø = 0

ø = x
ø' = y

Initial conditions on x and y 

x' = y 
y'  = - omega*x

So 

x' = 0*x 	+ 1* y
y' = -omega*x + 0*y

M = 0 		1
	-omega 	0
"""


"""
def harmonicOscillator(omega = 2):
	M = np.array([[0, 1], [-omega, 0]])
	return lambda t, y: M@y

t, yArray = rungeKuttaODE(f = harmonicOscillator(), y0 = np.array([1,0]), t0 = 0, tmax = 6, N = 1000)

plt.plot(yArray[:, 0], yArray[:, 1], color = 'green')
plt.show()


def dampedHarmonicOscillator(q, a, wD):
	return lambda t, y: np.array([y[1],-np.sin(y[0]) - (1/q)*(y[1]) + a*np.cos(y[2]), wD ])


t, yArray = rungeKuttaODE(f = dampedHarmonicOscillator(2, 0.9, 2/3), y0 = np.array([1,0,0]), t0 = 0, tmax = 500, N = 10000)	

plt.scatter(yArray[:, 0], yArray[:, 1], color = 'blue', s = 2)
plt.show()


"""



