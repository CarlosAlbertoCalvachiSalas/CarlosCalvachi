import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('gaussianDerivative.csv')


x = data['x']
yPrime = data['yPrime']


plt.plot(x, yPrime, color = 'blue')
plt.title(r'd/dx($\exp{(-x^2)}$)')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
plt.close()