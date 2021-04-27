import numpy as np
import matplotlib
import scipy
from matplotlib import pyplot as plt

x0 = -3
x1 = -2

def f(x):
    return (x**2)*np.sin(x)

def tolerance_test(x0,xprime):
    epsilon = 10**(-8)
    if abs(xprime-x0) < epsilon:
        return True
    else:
        return False

def newton(x0,x1):
    y0 = f(x0)
    y1 = f(x1)
    xprime = x1 - y1*((x1-x0)/(y1-y0))
    return xprime

def newton_method(x0,x1):
    xprime = x1
    j = 0
    while tolerance_test(x0,xprime) == False:
        xprime = newton(x0,x1)
        x0 = x1
        x1 = xprime
        j += 1
    return xprime, j

x = np.linspace(-5,5,100)
plt.plot(x,f(x))
plt.grid(True)
plt.show()

root, iterations = newton_method(x0,x1)
print(root)
print(str(iterations) + ' iterations')
print(np.pi)
