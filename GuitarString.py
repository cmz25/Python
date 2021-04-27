import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as animation

h = 2.0
l = 10.0   #you can format these to change the function
d = 5.0
N = 5

fig = plt.figure()
ax = plt.axes(xlim=(0, l), ylim=(-h-1, h+1))
line, = ax.plot([], [], color='brown', lw=0.75)
A = (h*(2**(1/2))*l**(3/2))/(2.25)

def C(n):
    return ((np.sin((n*np.pi*d)/l))/(n*np.pi)**2)*((1/d)+(1/(l-d)))+(np.sin((n*np.pi*d)/l)/(n*np.pi))*((1/(l-d))-(1/l)-(d/(l**2-(l*d))))

def init():
    line.set_data([], [])
    return line,

def animate(i):
    xpoints = np.linspace(0, l, 1000)
    y = np.zeros(len(xpoints),float)
    for n in range(1,N):
        j = 0
        for x in xpoints:
            y[j] += C(n)*A*np.sin((n*np.pi*x)/l)*np.cos((n*np.pi*0.2*i)/l)
            j+=1
    line.set_data(xpoints,y)
    return line,

anim = animation.FuncAnimation(fig,animate,init_func=init,frames=100,interval=20,blit=True)
plt.show()
