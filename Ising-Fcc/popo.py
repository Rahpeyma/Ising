import numpy as np
import matplotlib.pyplot as plt
import numba as nb
from timeit import default_timer as timer
# Define lattice parameters
a = 1.0  # lattice constant
N = 12  # number of unit cells per side
J1 = 1.0  # nearest-neighbor interaction strength
J2 = -0.1  # next-nearest-neighbor interaction strength
kB = 1.0  # Boltzmann constant
l=4


@nb.njit()
def Hamiltonian(lattice,i,j,k,o):
    #o=1
    if(o==1):
        energy=-J1*lattice[i,j,k,o]*(lattice[i,j,k,o-1]+lattice[i,j,k,o+1]+lattice[i,j,k,o+2]
                                  +lattice[(i+1)%N,j,k,o-1]+lattice[(i+1)%N,j,k,o+2]
                                  +lattice[i,(j+1)%N,k,o-1]+lattice[i,(j+1)%N,k,o+1]
                                  +lattice[i,j,(k-1)%N,o+1]+lattice[i,j,(k-1)%N,o+2]
                                  +lattice[(i+1)%N,(j+1)%N,k,o-1]
                                  +lattice[i,(j+1)%N,(k-1)%N,o+1]
                                  +lattice[(i+1)%N,j,(k-1)%N,o+2]
                                  )
    # o=2 
    elif(o==2):
        energy=-J1*lattice[i,j,k,o]*(lattice[i,j,k,o-2]+lattice[i,j,k,o-1]+lattice[i,j,k,o+1]
                                  +lattice[(i+1)%N,j,k,o-2]+lattice[(i+1)%N,j,k,o+1]
                                  +lattice[i,(j+1)%N,k,o-2]+lattice[i,(j+1)%N,k,o]
                                  +lattice[i,j,(k-1)%N,o]+lattice[i,j,(k-1)%N,o+1]
                                  +lattice[(i+1)%N,(j+1)%N,k,o-2]
                                  +lattice[i,(j+1)%N,(k-1)%N,o]
                                  +lattice[(i+1)%N,j,(k-1)%N,o+1]
                                  )
    # o=0 
    elif(o==0):
        energy=-J1*lattice[i,j,k,o]*(lattice[i,j,k,o+1]+lattice[i,j,k,o+2]+lattice[i,j,k,o+3]
                                  +lattice[(i-1)%N,j,k,o+1]+lattice[(i-1)%N,j,k,o+2]
                                  +lattice[i,(j-1)%N,k,o+1]+lattice[i,(j-1)%N,k,o+3]
                                  +lattice[i,j,(k-1)%N,o+2]+lattice[i,j,(k-1)%N,o+3]
                                  +lattice[(i-1)%N,(j-1)%N,k,o+1]
                                  +lattice[(i-1)%N,j,(k-1)%N,o+2]
                                  +lattice[i,(j-1)%N,(k-1)%N,o+3]
                                  )
    # o=3
    elif(o==3): 
        energy=-J1*lattice[i,j,k,o]*(lattice[i,j,k,o-3]+lattice[i,j,k,o-2]+lattice[i,j,k,o-1]
                                  +lattice[(i-1)%N,j,k,o-2]+lattice[(i-1)%N,j,k,o-1]
                                  +lattice[i,(j+1)%N,k,o-3]+lattice[i,(j+1)%N,k,o-1]
                                  +lattice[i,j,(k+1)%N,o-3]+lattice[i,j,(k-1)%N,o-2]
                                  +lattice[(i-1)%N,(j+1)%N,k,o-1]
                                  +lattice[(i-1)%N,j,(k+1)%N,o-2]
                                  +lattice[i,(j+1)%N,(k+1)%N,o-3]
                                  )
    energy+=-J2*lattice[i,j,k,o]*(lattice[(i+1)%N,j,k,o]+lattice[i,(j+1)%N,k,o]+lattice[i,j,(k+1)%N,o]
                                +lattice[(i-1)%N,j,k,o]+lattice[i,(j-1)%N,k,o]+lattice[i,j,(k-1)%N,o])
    return energy

    # Compute the energy change due to flip of spin i
    delta_E = 0
    for j in range(pos.shape[0]):
        if i != j:
            delta = pos[i] - pos[j]
            r = np.linalg.norm(delta)
            r = r - L*np.round(r/L) # Periodic boundary condition
            r=np.abs(r)
            mask_nn = r < a_1  # nearest-neighbor mask
            delta_E +=   (J1*spins[i]*spins[j]*mask_nn)
    
    # Flip the spin with probability exp(-dE/(kB*T))
    if np.random.rand() < np.exp(-2*delta_E / (kB * T)):
        spins[i] *= -1
    
    return spins

@nb.njit()
def metropolis_update(lattice, T):
    i = np.random.randint(0, N)
    j = np.random.randint(0, N)
    k = np.random.randint(0, N)
    o = np.random.randint(0, l)
    dE = Hamiltonian(lattice,i,j,k,o)
    if np.random.rand() < np.exp(2*dE/(kB*T)):
        lattice[i,j,k,o] *= -1

    return lattice

# Define function for Monte Carlo simulation
@nb.njit()
def monte_carlo(T, nsteps, lattice):

    for step in range(nsteps):
        lattice = metropolis_update(lattice, T)
        # print(spins)
    M = np.sum(lattice)
    # print(spins)
    result=M/(12*12*12*4)
    return result
    
    
@nb.njit()
def main():

    lattice=np.ones((N,N,N,l))
    # print(spins)
    T_range = np.linspace(0.1, 7, 10)
    res = np.zeros(len(T_range))
    for i in range(len(T_range)):
        r = monte_carlo(T_range[i], 10**(5), lattice)
        # print(spins)
        res[i] = r
    return res
    # plt.plot(T_range,res)
    # plt.show()


if __name__ == '__main__':
    T_range = np.linspace(0.1, 7, 10)
    start_time = timer()
    plt.plot(T_range,np.abs(main()))
    plt.show()
    end_time = timer()
    print("Elapsed time: %g seconds" % (end_time - start_time))
