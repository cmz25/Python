import numpy as np
import matplotlib
from matplotlib import pyplot as plt


L = 1
def phi1(i,N):
    x = (i/N)*L
    return x**2

def phi2(i,N):
    x = (i/N)*L
    return x**4-x**2+x

def Rho_func(y,x):
    return np.exp(-((x-0.45)**2+(y-0.65)**2)/0.02)

def charge_dist(N):
    rho = np.zeros((N+1,N+1),float)
    for i in range(N+1):
        for j in range(N+1):
            rho[i,j] = Rho_func(i*(L/N),j*(L/N))
    return rho

def boundary(N):
    potential = np.zeros((N+1,N+1),float)
    potential[1:N,0] = 0 #x=0
    potential[1:N,N] = 1 #x=L
    V2 = np.zeros((N-1),float)
    V3 = np.zeros((N-1),float)
    for i in range(1,N):
        V2[i-1] = phi1(i,N)
        V3[i-1] = phi2(i,N)
    potential[N,1:N] = V3[:] #y=0
    potential[0,1:N] = V2[:] #y=L
    return potential

def relaxation(P,N):
    P_prime = np.zeros((N+1,N+1),float)
    P_prime[:,:] = P[:,:]
    rho = charge_dist(N)
    dx = L/N
    for i in range(1,N):
        for j in range(1,N):
            P_prime[i,j] = (P[i-1,j]+P[i+1,j]+P[i,j-1]+P[i,j+1])/4 + (np.pi*dx**2)*rho[i,j]
    return P_prime

def equate(N,potential):
    P = boundary(2*N)
    for i in range(1,N):
        for j in range(1,N):
            P[2*i,2*j] = potential[i,j]
    return P

def interpolate(N,P):
    for i in range(2,N,2): #horizontal
        for j in range(1,N,2):
            P[i,j] = (P[i,j+1]+P[i,j-1])/2
    for i in range(1,N,2): #vertical
        for j in range(2,N,2):
            P[i,j] = (P[i+1,j]+P[i-1,j])/2
    for i in range(1,N,2): #stencil
        for j in range(1,N,2):
            P[i,j] = (P[i+1,j]+P[i-1,j]+P[i,j+1]+P[i,j-1])/4
    return P

def plot_potential(N,potential):
    plt.imshow(potential,cmap='Greys',origin='lower')
    cbar = plt.colorbar()
    cbar.set_label('Potential Φ')
    x = np.array([0,0.5,1])
    y = np.array([0,0.5,1])
    x_positions = np.array([0,N//2,N],float)
    y_positions = np.array([0,N//2,N],float)
    plt.xticks(x_positions, x)
    plt.yticks(y_positions, y)
    x_contour = np.arange(0,N+1,1)
    y_contour = np.arange(0,N+1,1)
    plt.contour(x_contour,y_contour,potential, levels = 20, vmin = 0, vmax = 1, cmap = 'YlOrRd')
    cbar = plt.colorbar()
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.title(f'Electric Potential Φ Contour N = {N} Multi-Grid Method')
    plt.show()

def multi_grid(dx):
    N = 2
    max = int(-np.log2(dx)-1)
    potential = boundary(N)
    potential_prime = relaxation(potential,N)
    plot_potential(N,potential)
    for i in range(max):
        potential_new = equate(N,potential_prime)
        print(f'\nThe potential at [0.5,0.5] for a discretization of N = {N} is {potential_prime[N//2,N//2]}')
        N = 2*N
        potential = interpolate(N,potential_new)
        potential_prime = relaxation(potential,N)
        plot_potential(N,potential)
    print(f'\nThe potential at [0.5,0.5] for a discretization of N = {N} is {potential_prime[N//2,N//2]}')

dx_final = 2**(-10)
multi_grid(dx_final)
