import numpy as np
import matplotlib.pyplot as plt 
import os.path as path
import wget

def LoadFiles(file, url):
    if not path.exists(file):
        Path_ = wget.download(url,file)
        print('File loaded')
    else:
        Path_ = file
        
    return Path_


LoadFiles('data.txt', 'https://raw.githubusercontent.com/asegura4488/Database/main/MetodosComputacionalesReforma/EnergiaPotencialGas2D.txt')

data = np.loadtxt('data.txt')

time   = data[:, 0]  
energy = data[:, 1]

normalizedEnergy = energy - np.average(energy)

fourierTransform = np.fft.fft(normalizedEnergy)
frequencies      = np.fft.fftfreq(len(time), d = time[1] - time[0])

modifiedFourierTransform = 2*(frequencies > 0)*fourierTransform

maxFourierAmplitude = np.max(modifiedFourierTransform)

mainFrequency = frequencies[modifiedFourierTransform == maxFourierAmplitude][0]
 
modifiedFourierTransform  = modifiedFourierTransform * (modifiedFourierTransform == maxFourierAmplitude)

mainOscillation = np.fft.ifft(modifiedFourierTransform)

plt.title('Tiempo libre {:.2f} pasos temporales'.format(1000/mainFrequency))
plt.plot(time*1000, energy, linestyle = '--', label = 'Serie Original')
plt.scatter(time*1000, np.real(mainOscillation) + np.average(energy), label = 'Fundamental', s = 7, color = 'darkorange')
plt.xlabel('t[s]')
plt.ylabel('E[J]')
plt.legend()
plt.show()
