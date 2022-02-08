import numpy as np 
import matplotlib.pyplot as plt

data = np.loadtxt('ManchasSolares.dat')

filteredData = data[data[:, 0] >= 1900]

years = filteredData[:, 0]

months = filteredData[:, 1]

time = years + months/12

spots = filteredData[:, 3]

fft = (np.fft.fft(spots - np.average(spots)))*2
freqfft = np.fft.fftfreq(len(time), d = 1/12)

positiveFFT 	= fft[freqfft > 0]
positiveFreqFFT = freqfft[freqfft > 0]

maximumFrequency = positiveFreqFFT[np.max(np.abs(positiveFFT)) == np.abs(positiveFFT)][0]
amplitude = positiveFFT[np.max(np.abs(positiveFFT)) == np.abs(positiveFFT)][0] 

amplitude = np.abs(amplitude)

fftFiltered = fft * (freqfft == maximumFrequency)

period = 1/maximumFrequency

principalOscillation = np.fft.ifft(fftFiltered)

plt.plot(time, spots, label = 'Datos')
plt.plot(time, np.real(principalOscillation) + np.average(spots), color = 'red', label = 'Frecuencia dominante')
plt.title('Periodo {:.2f} en años'.format(period))
plt.ylabel('Norma FFT')
plt.xlabel('Frencuencia[1/mes]')
plt.legend(loc = 'upper right')
plt.show()




"""
def getFFTPerYear(year):
	yearFilter = filteredData[filteredData[:, 0] == year]

	yearSpots = yearFilter[:, 2]

	normFTT = np.abs(np.fft.fft(yearSpots))
	freqFTT = np.fft.fftfreq(len(yearSpots))

	maxNorm = np.max(normFTT[freqFTT >= 0])
	maxFreq = freqFTT[normFTT == maxNorm][0]

	print(year)

	return maxFreq



def sumTotalSpots(year):
	yearFilter = filteredData[filteredData[:, 0] == year]

	yearSpots = yearFilter[:, 2]

	return np.sum(yearSpots)



sumTotalSpots = np.vectorize(sumTotalSpots)(np.unique(years))


plt.plot(np.unique(years), sumTotalSpots)
plt.show()



#normsFFT = np.vectorize(getFFTPerYear)(np.unique(years))




plt.plot(np.unique(years), normsFFT)
plt.show()
"""


"""

x = np.linspace(0, 10, 100000)
y = np.sin(2*x) + np.sin(3*x) + np.sin(100*x) + np.sin(200*x)

plt.plot(x, y, color = 'blue')
plt.show()
plt.close()

normFTT = np.abs(np.fft.fft(y))
freqFTT = np.fft.fftfreq(len(x), d = x[1] - x[0])



maxNorm = np.max(normFTT[freqFTT >= 0])
maxFreq = freqFTT[normFTT == maxNorm][0]


plt.scatter(freqFTT, normFTT, color = 'blue')
plt.show()
plt.close()

"""




