import numpy as np
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm 
import matplotlib.animation as anim

sys.setrecursionlimit(100000)

# Euler Rule

# Equations of the form y'(t) = f(t) given initial condition 
# y(0) = a, with a as a constant.


def rungeKuttaODE2(f, y0, N):
	t = np.arange(N+1)
	h = t[-1]/(N)
	y = [y0]

	for x in tqdm(range(1, N+1)):
		yn_1 = y[x-1]
		k1 = h*f(t[x-1], yn_1)
		k2 = h*f(t[x-1]+(h/2), yn_1 + (k1/2))
		k3 = h*f(t[x-1]+(h/2), yn_1 + (k2/2))
		k4 = h*f(t[x-1]+h, yn_1 + k3)

		yn = yn_1 + ((k1+2*k2+2*k3+k4)/6)

		y.append(yn)

	y = np.array(y)

	return t, y

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


G = 6.67E-11
mT = 5.9736E24
rT = 6.3781E6
mL = 0.07349E24
rL = 1.7374E6
d  = 3.844E8
omega = 2.6617E-6

Delta = (G*mT)/(d**3)
mu = mL/mT

def lunarRocketTrayectoryEquation(t, y):
	r 		= y[0]
	phi 	= y[1]
	pR 		= y[2]
	pPhi 	= y[3]

	rprime = np.sqrt(1 + r**2 - 2*r*np.cos(phi - omega*t))
		
	return np.array([pR, pPhi/(r**2) , (pPhi**2)/(r**3) - Delta*((1/r**2) + ((mu/(rprime**3))*(r - np.cos(phi - omega*t)))),  - ((Delta*mu*r)/(rprime**3))*np.sin(phi - omega*t)])	

def lunarRocketTrayetory(nu, phi0, theta):
	v0 = np.sqrt((nu*G*mT)/rT)/d

	r0 = rT/d
	phi0 = phi0 
	pr0 = v0*np.cos(theta - phi0)
	pphi0 = r0*v0*np.sin(theta - phi0)

	t, yArray = rungeKuttaODE2(f = lunarRocketTrayectoryEquation, y0 = np.array([r0,phi0,pr0,pphi0]), N = 3600*24*2 )	

	#t, yArray = rungeKuttaODE(f = lunarRocketTrayectoryEquation, y0 = np.array([r0,phi0,pr0,pphi0]), t0 = 0, tmax = 3600*24*2, N = 10000)	

	x = yArray[:, 0]*np.cos(yArray[:, 1])
	y = yArray[:, 0]*np.sin(yArray[:, 1])

	xMoon = np.cos(omega*t)
	yMoon = np.sin(omega*t)

	return t, x, y, xMoon, yMoon

t, x, y, xMoon, yMoon = lunarRocketTrayetory(nu = 2.0035, phi0 = 0, theta = 0.0872*np.pi)


t = t[::1000]
x = x[::1000]
y = y[::1000]
xMoon = xMoon[::1000]
yMoon = yMoon[::1000]

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)

def init():
    ax.set_xlabel('x/d')
    ax.set_ylabel('y/d')

def Update(i):
    
    plot = ax.clear()
    init()
    ax.set_title(r'$t=%.2f \ h$ LunarRocket ' %(t[i] / 3600), fontsize=15)
    
    circleEarth = plt.Circle((0,0), rT/d, color='blue', fill=True, label = 'Tierra')
    circleMoon  = plt.Circle((xMoon[i],yMoon[i]), rL/d, color='grey', fill=True, label = 'Luna')
    circleSpacecraft = plt.Circle((x[i],y[i]), (rL/d)/10, color='red', fill=True, label = 'Nave')

    ax.add_patch(circleEarth)
    ax.add_patch(circleMoon)
    ax.add_patch(circleSpacecraft)

 

    plt.legend(loc = 'upper right')
        
    return plot

Animation = anim.FuncAnimation(fig, Update, frames=len(t), init_func=init, interval = 10, repeat = False)

plt.show()


