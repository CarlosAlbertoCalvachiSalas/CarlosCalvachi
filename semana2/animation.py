import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera

t = np.linspace(-np.pi, np.pi, 10000)

def seriesTerm(n, t):
	return 2 * ((((-1)**(n-1)) * np.sin(n*t))/n)

def seriesTermVect(n,t):
	return np.vectorize(seriesTerm, excluded = ['t'], signature = '()->(n)', cache = True)(n = n, t = t)

harmonicNums = np.arange(50) + 1
y = np.cumsum(seriesTermVect(n = harmonicNums , t = t), axis = 0)

fig = plt.figure()
plt.xlim(-np.pi, np.pi)
plt.title('Armónicos de la función ' + r'$f(t) = t$')

camera = Camera(fig)

for x in range(len(y)):
	plt.plot(t, y[x], color = 'blue')
	camera.snap()

animation = camera.animate(interval = 200)
plt.show()
