import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from tqdm import tqdm 

k = 5 # N/m
l = 3 #
m = 2 # kg
g = 9.8 # m/s2
r0 = 15 #m
theta0 = np.pi/8
r0prime = 0
theta0prime = 0



def pendular(beta = 0):

	N = 2001

	t = np.linspace(0, 30, N)

	h = t[1] - t[0]

	# Hay dos opciones. O se inicia con las posiciones x y y 
	# o se inicia con las variables r y theta 
	# Primer paso puede ser Euler

	r = np.zeros(N)
	theta = np.zeros(N)
	dotr = np.zeros(N)
	dottheta = np.zeros(N)
	
	# Se disponen las aceleraciones

	def ddotr(r, theta, rprime, thetaprime, beta = beta):
		return (k/m)*(l-r) + g*np.cos(theta) + r*(thetaprime**2) - beta*rprime

	def ddottheta(r, theta, rprime, thetaprime, beta = beta):
		return (-g/r)*np.sin(theta) - (2/r)*(rprime*thetaprime) - beta*thetaprime

	# Se inicializan las condiciones inicales

	r[0] = r0
	theta[0] = theta0
	dotr[0] = r0prime
	dottheta[0] = theta0prime

	# Se utiliza el método de Euler para el siguiente paso 

	dotr[1] 	= dotr[0] + h*ddotr(r[0], theta[0], dotr[0], dottheta[0])
	dottheta[1] = dottheta[0] + h*ddottheta(r[0], theta[0], dotr[0], dottheta[0])
	r[1]		= r[0] + h*dotr[0]
	theta[1]    = theta[0] + h*dottheta[0] 

	# Se corre el inicializador 

	for x in range(2, N):

		arx = ddotr(r[x-1], theta[x-1], dotr[x-1], dottheta[x-1])
		arx_1 = ddotr(r[x-2], theta[x-2], dotr[x-2], dottheta[x-2])

		athetax = ddottheta(r[x-1], theta[x-1], dotr[x-1], dottheta[x-1])
		athetax_1 = ddottheta(r[x-2], theta[x-2], dotr[x-2], dottheta[x-2])

		r[x] = r[x-1] + h*dotr[x-1] + (1/6)*(4*arx - arx_1)*(h**2)
		dotr[x] = dotr[x-1] + (1/2)*(3*arx - arx_1)*h

		theta[x] = theta[x-1] + h*dottheta[x-1] + (1/6)*(4*athetax - athetax_1)*(h**2)
		dottheta[x] = dottheta[x-1] + (1/2)*(3*athetax - athetax_1)*h

	# Se corre el corrector

	for x in range(2, N):

		arx1 = ddotr(r[x], theta[x], dotr[x], dottheta[x])
		arx = ddotr(r[x-1], theta[x-1], dotr[x-1], dottheta[x-1])
		arx_1 = ddotr(r[x-2], theta[x-2], dotr[x-2], dottheta[x-2])

		athetax1 = ddottheta(r[x], theta[x], dotr[x], dottheta[x])
		athetax = ddottheta(r[x-1], theta[x-1], dotr[x-1], dottheta[x-1])
		athetax_1 = ddottheta(r[x-2], theta[x-2], dotr[x-2], dottheta[x-2])

		dotr[x] = dotr[x-1] + (1/12)*(5*arx1 + 8*arx - arx_1)*h
		dottheta[x] = dottheta[x-1] + (1/12)*(5*athetax1 + 8*athetax - athetax_1)*h

	# Se produce una reducción de factor de 20

	r = r[::20]
	theta = theta[::20]
	t = t[::20]
	# Se calculan las posiciones x, y

	x = r*np.cos(theta - np.pi/2)
	y = r*np.sin(theta - np.pi/2)


	fig = plt.figure(figsize=(15,10))
	ax1 = fig.add_subplot(1,2,1, projection='polar')
	ax2 = fig.add_subplot(1,2,2,)

	def init():
		ax2.set_xlim(-20,20)
		ax2.set_ylim(-20,20)

	def Update(i):
		ax2.clear()
		ax2.grid()
		ax2.axis('equal')
		init()
		ax2.set_title(r't=%.2f seconds ' %t[i] + r'$\beta = $' + str(beta), fontsize=15)
		circleBall = plt.Circle((x[i], y[i]), 2, color='black', fill=True, label = 'Tierra')

		ax2.add_patch(circleBall)
		ax2.arrow( 0., 0., x[i], y[i], head_width=0.5, color = 'r')
		ax1.scatter(theta[i], r[i])

	Animation = anim.FuncAnimation(fig,Update, frames=len(t), init_func=init, interval = 100, repeat= False)

	plt.show()

print('Oscilador pendular sin fricción')
pendular(beta = 0)
print('Oscilador pendular con fricción')
pendular(beta = 0.2)
