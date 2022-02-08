import numpy as np 
import matplotlib.pyplot as plt

def function(x):
	return np.exp(-x)*np.sin(x)

def exactDerivative(x):
	return np.exp(-x)*np.cos(x) - np.exp(-x)*np.sin(x)

def centralDerivative(x, h = 0.5):
	return (function(x + h) - function(x - h))/(2*h)

def spectralDerivative(signal, h = 0.5):
	fourierTransform = np.fft.fft(signal)
	frecuencies = np.fft.fftfreq(len(signal), d = h)
	derivative 	= fourierTransform*(1j)*2*frecuencies* 2 * np.pi
	derivative = derivative * (frecuencies >= 0)
	spectralDerivative = np.fft.ifft(derivative)

	return np.real(spectralDerivative)

x = np.arange(0, 2*np.pi, 0.5)

exactDerivative = np.vectorize(exactDerivative)(x)
centralDerivative = np.vectorize(centralDerivative)(x)

signal = np.vectorize(function)(x)

plt.plot(x, signal, color = 'blue', label = 'función')
plt.plot(x, exactDerivative, color = 'green', label = 'exacta')
plt.plot(x, centralDerivative, color = 'brown', label = 'central')
plt.plot(x, spectralDerivative(signal), color = 'red', label = 'espectral')
plt.title(r'Derivadas de la función $e^{-x}\sin{x}$ con $\Delta x = 0.5$')
plt.legend()
plt.show()

print('\n Si se comparan las expresiones, la derivada exacta está mas cerca de la derivada central. Por contrastre, la derivada espectral es más fluctuante, especialmente al final del intérvalo. Esto se debe probablemente a que se está aplicando la transformada inversa de un espectro de frecuencias que, en el dominio del tiempo, se puede ver como una suma de varias funciones oscillantes, por eso la forma aproximada de zig zag')

def function(x):
	return np.exp(-x)*np.sin(x)

def exactDerivative(x):
	return np.exp(-x)*np.cos(x) - np.exp(-x)*np.sin(x)

def centralDerivative(x, h = 0.05):
	return (function(x + h) - function(x - h))/(2*h)

def spectralDerivative(signal, h = 0.05):
	fourierTransform = np.fft.fft(signal)
	frecuencies = np.fft.fftfreq(len(signal), d = h)
	derivative 	= fourierTransform*(1j)*2*frecuencies* 2 * np.pi
	derivative = derivative * (frecuencies >= 0)
	spectralDerivative = np.fft.ifft(derivative)

	return np.real(spectralDerivative)

x = np.arange(0, 2*np.pi, 0.05)

exactDerivative = np.vectorize(exactDerivative)(x)
centralDerivative = np.vectorize(centralDerivative)(x)

signal = np.vectorize(function)(x)

plt.plot(x, signal, color = 'blue', label = 'función')
plt.plot(x, exactDerivative, color = 'green', label = 'exacta')
plt.plot(x, centralDerivative, color = 'brown', label = 'central')
plt.plot(x, spectralDerivative(signal), color = 'red', label = 'espectral')
plt.title(r'Derivadas de la función $e^{-x}\sin{x}$ con $\Delta x = 0.05$')
plt.legend()
plt.show()

print(' \n Si se aumenta la resolución, se observa que la derivada espectral aunque fluctua por las razones presentadas, converge a la derivada real. Esto es un test de consistencia')


