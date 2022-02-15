import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('gaussianDerivative.txt')

x = data[:, 0]
yPrime = data[:, 1]


plt.plot(x, yPrime, color = 'blue')
plt.title(r'd/dx($\exp{(-x^2)}$)')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
plt.close()