#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time

#--- constants ---
# units: eV, A, S
HBAR = 6.582119569E-16 #eV*s
REST_ENERGY_E = 0.51999895069E6 #eV
C = 2.997924580E18 # A/S
M_E = REST_ENERGY_E/(C**2)

#--- parameters ---
#size of simulation box
l = 300  #A
N = 2000 #number of points in grid
delta_x = l/(N-1)

#--- potentials ---
def qho_potential(x) -> float:
    hbar_omega = 10E-3
    return (1/2)*((hbar_omega/HBAR)**2)*M_E*(x**2)

def finite_well_potential(x) -> float:
    V_height = 5
    well_length = 10

    if x < well_length/2 and x > -well_length/2: 
        return 0
    else: 
        return V_height

# --- functions ---
def test_normalization(psi) -> bool: 

    # do the normalization integral
    int = 0
    for i in psi:
        int = int + (i**2)*delta_x

    # check if the integral -1 is close to zero (10^-18)
    if (int - 1) < 1E-18:
        return True
    else:
        # if the wave function wasnt normalized print error mesage
        print("the wave function is not normalized!")
        print("normalization integral: ", i)
        return False 
    
def solve_tise(l, N, potential, eigenstates=5, five_point = False): 

    if eigenstates >= N: 
        #print("warning: surpassed max limit of eigenstates (N-1)")
        #print(f"calculating {N-1} states instead")
        eigenstates = N - 1

    start_time = time.perf_counter()
    #discrete space
    x  = np.linspace(-l/2,l/2,N)
    dx = x[1] - x[0]
    
    #define kinetic term
    T = (HBAR**2)/(2*M_E*(dx**2))

    #pick potential and map it to an N sized array
    v = np.ones(N)
    for i in range(0,N):
        v[i] = potential(x[i])

    # build five point stencil hamiltonian (if needed)
    if five_point == True: 

        T_diag2 = (1/12) * T
        T_diag1 = (-4/3) * T
        main_diag = v + (5/2)*T

        diagonals = [T_diag2, T_diag1, main_diag, T_diag1, T_diag2]
        offsets   = [-2, -1, 0, 1, 2]
        H = sp.sparse.diags(diagonals, offsets, shape=(N, N), format='csc')

        #call solver
        energies, psi = sp.sparse.linalg.eigsh(H, k=eigenstates, which='SA')
    else: 
        #build three point stencil hamiltonian
        D  = np.ones(N)
        SD = np.ones(N-1)

        # build hamiltonian 
        D  = v + 2*T
        SD = SD * -T

        #call solver
        energies, psi = sp.linalg.eigh_tridiagonal(D,SD, select = 'i', select_range = (0,eigenstates))

    #finish counting time
    end_time = time.perf_counter()
    time_taken = end_time - start_time

    return time_taken, energies, psi


# part 1)
# calculate ground state for different N 

time_vs_n = np.zeros((2,100)) 
print(" time[S] vs N vs Ground Energy[eV]")
for i,n in enumerate(range(10,5000,50)):
    print(i,n)
    time_vs_n[0,i], energy, wavefunctions = solve_tise(200,n,qho_potential,eigenstates=5,five_point=False)
    time_vs_n[1,i] = n

#%%
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(time_vs_n[1,:], time_vs_n[0,:], linewidth=2, label="time vs N")
ax.set_title('3 point stencil QHO')
ax.set_xlabel('N points in grid')
ax.set_ylabel('Time spent in calculations')
ax.set_ylim((0,0.012))
ax.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('t_vs_n_3point.png')
plt.show()
# notes and data
# Quantum Harmonic oscillator

# 3 point stencil time vs N (and energy) for the first 5 eigenstates
# time[S]  |   N   | Ground Energy[eV]
#----------------------------------------
# 0.14188     10     0.02091
# 0.00020     20     0.00593
# 0.00014     30     0.00462
# 0.00009     40     0.00471
# 0.00009     50     0.00482
# 0.00015     100    0.00496
# 0.00023     150    0.00498
# 0.00027     200    0.00499
# 0.00042     300    0.00500
# 0.00052     400    0.00500
# 0.00077     500    0.00500
# 0.00155     1000   0.00500
# 0.00272     2000   0.00500
# 0.00526     4000   0.00500
# 0.00716     6000   0.00500

# 5 point stencil time vs N (and energy) for the first 5 eigenstates
# time[S]  |   N   | Ground Energy[eV]
#----------------------------------------
# 0.15287     10     0.02096
# 0.00014     20     0.00613
# 0.00009     30     0.00491
# 0.00017     40     0.00493
# 0.00016     50     0.00497
# 0.00033     100    0.00500
# 0.00052     150    0.00500
# 0.00076     200    0.00500
# 0.00185     300    0.00500
# 0.00345     400    0.00500
# 0.01166     500    0.00500
# 0.13051     1000   0.00500
# 0.83912     2000   0.00500
# 9.08061     4000   0.00500
# 35.67184    6000   0.00500

# Finite square well
# 3 point stencil
# time[S]  |   N   | Ground Energy[eV]
#----------------------------------------
# 0.08620     10     0.18262
# 0.00013     20     0.28989
# 0.00009     30     0.20877
# 0.00009     40     0.24396
# 0.00011     50     0.26806
# 0.00017     100    0.27968
# 0.00026     150    0.25792
# 0.00031     200    0.26498
# 0.00045     300    0.27217
# 0.00058     400    0.26653
# 0.00071     500    0.27047
# 0.00141     1000   0.26733
# 0.00274     2000   0.26757
# 0.00540     4000   0.26769
# 0.00799     6000   0.26773

#five point stencil
# time[S]  |   N   | Ground Energy[eV]
#----------------------------------------
# 0.08826     10     0.21046
# 0.00012     20     0.30938
# 0.00009     30     0.21445
# 0.00013     40     0.24828
# 0.00012     50     0.27137
# 0.00030     100    0.28061
# 0.00053     150    0.25829
# 0.00080     200    0.26519
# 0.00198     300    0.27227
# 0.00368     400    0.26659
# 0.00763     500    0.27051
# 0.12142     1000   0.26734
# 0.90401     2000   0.26758
# 9.68930     4000   0.26769
# 35.78233    6000   0.26773


# plot_energy_levels(x,v,energies,indices=[0,1,2,3,4,5,6,7,8,9,10])
# plot_wavefunctions(x,psi,indices=[0,1,2,3,4,5],density=True)


# %%
