import numpy as np
import matplotlib
import scipy
from matplotlib import pyplot as plt
from numpy.fft import fft

N = 30
L = 1
h = L/N

def f(x):
    return (x-0.5)**2

def DFT():
    gamma_k = np.zeros(N,complex)
    for k in range(N):
        for n in range(N):
            gamma_k[k] += (1/N)*f(n*h)*np.exp((-2j*np.pi*k*n)/N)
    return gamma_k

def IDFT(gamma_k):
    y_n = np.zeros(N)
    for n in range(N):
        for k in range(N):
            y_n[n] += gamma_k[k]*np.exp((2j*np.pi*k*n)/N)
    return y_n

def DFT2():
    c = np.zeros(N//2+1,complex)
    for k in range(N//2+1):
        for n in range(N):
            c[k] += f(n*h)*(np.exp((-2j*np.pi*k*n)/N))
    return c

def IDFT2(c):
    y = np.zeros(N)
    for n in range(N):
        for k in range(N//2+1):
            y[n] += (1/N)*c[k]*((np.exp((2j*np.pi*k*n)/N))+(np.exp((2j*np.pi*k*n)/N)))
    return y

def convolution(c1,c2):
    c = np.zeros(N,complex)
    for k in range(N):
        c[k] = c1[k]*c2[k]
    return c

x = np.linspace(0,L,N)

gamma_k = DFT()
y_n = IDFT(gamma_k)
plt.plot(x,y_n,color = 'green')
plt.show()
