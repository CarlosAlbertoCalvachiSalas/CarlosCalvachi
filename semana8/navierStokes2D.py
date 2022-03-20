import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm 

N = 50
u = np.zeros([N + 1, N + 1])
w = np.zeros([N + 1, N + 1])
h = 0.02
V0 = 1
nu = 0.2
R = (V0*h)/nu

omega = 0.5
NInstances = 1000

#	Obstacule

x0 = 5
y0 = 5

x1 = 25
y1 = 30

u[x0:x1+1, y0:y1+1] = 0

w[x0, y0:y1+1] = (-2/(h**2))*(u[x0 - 1 , y0:y1+1] - u[x0, y0:y1+1])
w[x1, y0:y1+1] = (-2/(h**2))*(u[x1 + 1 , y0:y1+1] - u[x1, y0:y1+1])
w[x0:x1+1, y1] = (-2/(h**2))*(u[x0:x1+1, y1+1] - u[x0:x1+1, y1])
w[x0:x1+1, y0] = (-2/(h**2))*(u[x0:x1+1, y0-1] - u[x0:x1+1, y0])

w[x0+1:x1, y0+1:y1] = 0


def relax():
	#beam()
	for x in range(1, N):
		for y in range(1, N):
			nextu = (1/4)*(u[x+1, y] + u[x-1, y] + u[x, y+1] + u[x, y-1] + (h**2)*w[x,y])
			r1 = omega*(nextu - u[x,y])
			u[x, y] = u[x, y] + r1

	for x in range(1, N):
		for y in range(1, N):
			homow = (1/4)*(w[x+1, y] + w[x-1, y] + w[x, y+1] + w[x, y-1])
			cross1 = (u[x, y + 1] - u[x, y-1])*(w[x+1, y] - w[x-1, y])
			cross2 =(u[x+1, y] - u[x-1, y])*(w[x, y + 1] - w[x, y-1])
			newW = homow + (R/16)*(cross2 - cross1)

			r2 = omega*(newW -  w[x,y])
			w[x, y] = w[x, y] + r2

	# Obstacule 
	u[x0:x1+1, y0:y1+1] = 0

	# Inlet
	u[1] = u[0]

	# Surface
	u[:, N] = u[:, N-1] + h*V0

	# Outlet

	u[N] = u[N-1]
	w[N] = w[N-1]

	# Beam

	w[x0, y0:y1+1] = (-2/(h**2))*(u[x0 - 1 , y0:y1+1] - u[x0, y0:y1+1])
	w[x1, y0:y1+1] = (-2/(h**2))*(u[x1 + 1 , y0:y1+1] - u[x1, y0:y1+1])
	w[x0:x1+1, y0] = (-2/(h**2))*(u[x0:x1+1, y0-1] - u[x0:x1+1, y0])
	w[x0:x1+1, y1] = (-2/(h**2))*(u[x0:x1+1, y1+1] - u[x0:x1+1, y1])
 	
	w[x0+1:x1, y0+1:y1] = 0
	
	return u, w 	

for x in tqdm(range(NInstances)):
	u, w = relax()

x = np.arange(0, N+1)
y = np.arange(0, N+1)

X, Y = np.meshgrid(x, y)

Z = u[X, Y]
plt.contour(X/50, Y/50, Z, levels = 30)
plt.colorbar()
plt.xlim(0,1)
plt.ylim(0,1)
plt.show()

Z = w[X, Y]
plt.contourf(X/50, Y/50, Z, levels = 20)
plt.colorbar()
plt.xlim(0,1)
plt.ylim(0,1)
plt.show()

left = u[:N-1, :] 
right = u[1:N, :]

midUx = (right - left) / (2*h)
uxDown = (u[1] - u[0])/h 
uxUp = (u[N] - u[N-1])/h 

ux = np.vstack((uxDown, midUx, uxUp))

down = u[:, :N-1]
up = u[:, 1:N]

midUy = (up - down) / (2*h)
uyDown = (u[:,1] - u[:,0])/h 
uyUp = (u[:,N] - u[:,N-1])/h 

uy = np.hstack((uyDown.reshape(-1,1), midUy, uyUp.reshape(-1,1)))

plt.streamplot(X/50, Y/50, uy[X, Y], -ux[X, Y], color = 'blue', density = 2)
plt.title('Velocity Field')
plt.xlim(0,1)
plt.ylim(0,1)
plt.show()

