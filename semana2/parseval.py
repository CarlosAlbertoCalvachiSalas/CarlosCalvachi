import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


print('Se utiliza la relación de parseval, y en particular, la igualdad de la integral de la ecuación (26) de la solución de la parte teórica')
def funcionParseval(x):
	return (1/np.pi)*(((1/12)*(x**3 - (((np.pi)**2) * x) ))**2)

integral, error = integrate.quad(funcionParseval, -np.pi, np.pi)

print('Integral numérica: ', integral, 'error: ', error)

print('Valor exacto', np.pi**6/945)

print('Error numérico contra exacto', np.abs(integral - np.pi**6/945))

t = np.linspace(-np.pi, np.pi, 10000)
signal = t**2 - np.average(t**2)

frecuencies = np.fft.fftfreq(len(signal), d = t[1] - t[0])

fourierTransform = (np.fft.fft(signal) * 2) 

fourierTransform[frecuencies > 0] = (1/(2*np.pi*frecuencies[frecuencies > 0]*(1j)))*fourierTransform[frecuencies > 0]

fourierTransform = fourierTransform * (frecuencies > 0) 

bn = (np.sum(np.real(fourierTransform))/16) 


plt.plot(t, np.real(np.fft.ifft(fourierTransform)) )
#plt.plot(t, t**3/3 -np.average(t**2)*t , color = 'red')
plt.show()