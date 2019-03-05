import numpy as np
import matplotlib.pyplot as plt
import scipy.special
import scipy.interpolate
import scipy.integrate
import copy

# Import the input file (Input file assumed to be in current working directory)
#import sys
#import os
#workdir = os.getcwd()
#sys.path.append('{}'.format(workdir))

###############################################################################
#                                                                             #
# Compute the analytical analytical solution of 1-D radiative transfer in an  #
# isotropiclly scattering media with black or diffuse gray walls.             #
#                                                                             #
#-----------------------------------------------------------------------------#
#                                                                             #
# Based on the solution method presented in chapter 14 section 4 of:          #
# Radiative Heat Transfer 3rd. Edition,                                       #
# Michal F. Modest,                                                           #
# Academic Press, 2013                                                        #
#                                                                             #
###############################################################################

##############################################################################
#                                Constants                                   #
##############################################################################

# Stefan Boltzmann Constant - [W/(m^2*K^4)] 
sigma =  5.67036713*np.power(10.0,-8.0)

##############################################################################
#                                User Input                                  #
##############################################################################

# Size of domain
X_max = 1.0 # [m]

# Number of steps between 0 and X_max
Nx_analytical = 1001 # [points] 

# Wall temperatures
#T_left_wall  = 1.0 # [K]
#T_right_wall = 1.0 # [K]

#T_left_wall   = np.power(1/sigma, (1.0/4.0)) # [K]
T_left_wall   = 64.8 # [K]
#T_right_wall  = np.power(1/sigma, (1.0/4.0)) # [K]
T_right_wall  = 0.0 # [K]

# Wall Emsivities
epsilon_left_wall  = 1.0 # [unitless]
epsilon_right_wall = 1.0 # [unitless]

# Temperature of participating media - [K]
def T_gas(x):
    T_max = 0.0 # [K]
    #T_max = np.power(1/sigma, (1.0/4.0)) # [K]

    # Constant
    T =  T_max*np.ones([len(x)])

    ## Linear
    #T =  T_max + -1*((T_max-0)/X)*x

    return T # [K]

# Absorbivitiy of participating media - [1/m]
def sigma_a_gas(x):
    a_max = 1.0 # [1/m]

    # Constant
    sigma_a =  a_max*np.ones([len(x)]) 

    ## Linear
    #sigma_a = a_max + -1*((a_max - 0)/X)*x 

    ## Quadratic
    #sigma_a = (a_max/(X*X))*np.power(x-X,2)

    return sigma_a # [1/m]


# Scattering coefficent of participating media - [1/m]
def sigma_s_gas(x):
     sigma_a_max = 0.0  # [1/m]

     # Constant
     sigma_a = sigma_a_max *np.ones([len(x)])

     return sigma_a # [1/m]

# Name of output file
output_file_name = '1D_analytical_solution'     



##############################################################################
#                 Use user input to finish setting up the probem             #
##############################################################################

# Discritize spatial domain (x axis)
x = np.linspace(0.0, X_max, Nx_analytical)

# Discritize the input participating media properties on x
T = T_gas(x)                   # [K]
sigma_a = sigma_a_gas(x)       # [1/m]
sigma_s = sigma_s_gas(x)       # [1/m]
sigma_e = sigma_a + sigma_s    # [1/m]
omega = np.divide(sigma_s,sigma_e)   # [unitless]

# Optical Depth Units (Normal to Wall)
tau = np.zeros([len(x)])       # [unitless]
for i in range(0,len(x)):
    tau[i] = np.trapz(sigma_e[0:i+1],x[0:i+1])

# Print a summary of the input to the terminal
print('T_left_wall = {}'.format(T_left_wall))
print('T_right_wall = {}'.format(T_right_wall))
print('espilon_left_wall = {}'.format(epsilon_left_wall))
print('espilon_right_wall = {}'.format(epsilon_right_wall))
print('T = {}'.format(T))
print('sigma_a = {}'.format(sigma_a))
print('sigma_s = {}'.format(sigma_s))
print('sigma_e = {}'.format(sigma_e))
print('omega = {}'.format(omega))
print('tau = {}'.format(tau))

###############################################################################
##                              Assert the input is correct                   #
###############################################################################
#
## Make sure the walls are black 
#if MCRT_input.epsilon_left_wall(0,0) != 1:
#   print('ERROR')
#   print('MCRT_input.epsilon_left_wall(0,0) != 1')
#   exit()
#if MCRT_input.epsilon_right_wall(0,0) != 1:
#   print('ERROR')
#   print('MCRT_input.epsilon_right_wall(0,0) != 1')
#   exit()

##############################################################################
#                                 Solution                                   #
##############################################################################

# Black Body intensity of the gas
I_b = sigma*np.power(T,4)/(np.pi)

# Determine if the gas is scattering and take inital guess for the source term
if sum(sigma_s) == 0:
   scattering = False
   S = I_b
else:
   scattering = True
   S = np.multiply((1-omega),I_b)

# Radiocities of the walls
J = np.empty(2)
if (epsilon_left_wall == 1.0) and (epsilon_right_wall == 1.0): 
  # Both walls black
  J[0] = sigma*np.power(T_left_wall,4)   # Left  wall
  J[1] = sigma*np.power(T_right_wall,4)  # Right wall
elif (epsilon_left_wall == 1.0) and (epsilon_right_wall != 1.0) :
  # Left wall black right wall diffuse gray
  J[0] = sigma*np.power(T_left_wall,4)   # Left  wall
  J[1] = ( epsilon_right_wall/(1-epsilon_right_wall)*sigma*np.power(T_right_wall,4.0) + 2*J[0]*scipy.special.expn(3,tau[-1]) + 2*np.pi*np.trapz(S*scipy.special.expn(2, tau[-1]-tau) , tau) ) / (epsilon_right_wall/(1-epsilon_right_wall) + 2*scipy.special.expn(3,0.0) ) # Right Wall
elif (epsilon_left_wall != 1.0) and (epsilon_right_wall != 1.0) :
  # Left wall difuse gray right wall black
  J[1] = sigma*np.power(T_right_wall,4)  # Right wall
  J[0] = (epsilon_left_wall/(1-epsilon_left_wall)*sigma*np.power(T_left_wall,4.0) + 2*scipy.special.expn(3,tau[-1]) + 2*np.pi*np.trapz(S*scipy.special.expn(2, tau) , tau) ) / ( epsilon_left_wall/(1-epsilon_left_wall) + 2*scipy.special.expn(3,0.0) )  # Left  wall
elif (epsilon_left_wall != 1.0) and (epsilon_right_wall != 1.0) :
  # Left wall and right wall difuse gray
  A = np.array([[2*scipy.special.expn(3,0.0) + epsilon_left_wall/(1-epsilon_left_wall), -2*scipy.special.expn(3,tau[-1])], \
                [2*scipy.special.expn(3,tau[-1]),                                       -2*scipy.special.expn(3,0) - epsilon_right_wall/(1-epsilon_right_wall)] ])
  
  b = np.array([[epsilon_left_wall/(1-epsilon_left_wall)*sigma*np.power(T_left_wall,4) + 2*np.pi*np.trapz(S*scipy.special.expn(2, tau) , tau)], \
                [-1*epsilon_right_wall/(1-epsilon_right_wall)*sigma*np.power(T_right_wall,4) - 2*np.pi*np.trapz(S*scipy.special.expn(2, tau[-1]-tau) , tau)] ])
  
  J = np.linalg.solve(A,b)

# Arrays for the solution
G     = np.zeros([Nx_analytical])
G_old = np.zeros([Nx_analytical])
q     = np.zeros([Nx_analytical])

res = 1.0
while res > np.power(10.0,-6) :

    for i in range(0,len(tau)):

#        #######################################################################
#        # Evaluate the integrals in the G equation using the quad function
#        # Create interpolation function to use in the integral
#        S_interp = scipy.interpolate.interp1d(tau, S, kind='cubic')
#
#        # Check if the first point
#        if tau[i] == 0:
#            (integral_1, integral_1_error) = (0, 0) # Bounds make integral 0 in this case
#        else:
#            (integral_1, integral_1_error) = scipy.integrate.quad(lambda tau_p: S_interp(tau_p)*scipy.special.expn(1, tau[i]-tau_p),      0, tau[i])
#
#        # Check if the last point
#        if tau[i] == tau[-1]:
#            (integral_2, integral_2_error) = (0, 0) # Bounds make integral 0 in this case
#        else:
#            (integral_2, integral_2_error) = scipy.integrate.quad(lambda tau_p: S_interp(tau_p)*scipy.special.expn(1, tau_p-tau[i]), tau[i], tau[-1])
#
##        ## DEBUG ##
##        print(integral_1, integral_1_error)
##        print(integral_2, integral_2_error)
##        print('')
##        ## DEBUG ##
#
#        G[i] = 2*J[0]*scipy.special.expn(2,tau[i]) + \
#               2*J[1]*scipy.special.expn(2,tau[-1]-tau[i] ) + \
#               2*np.pi*integral_1 + \
#               2*np.pi*integral_2
#
#        #######################################################################

        #######################################################################
        # Approx G integrals by [0:i] and [i+1:] instad of [0:i+1] and [i:] to keep E1 from blowing up
        G[i] = 2*J[0]*scipy.special.expn(2,tau[i]) + \
               2*J[1]*scipy.special.expn(2,tau[-1]-tau[i] ) + \
               2*np.pi*np.trapz(S[0:i ]*scipy.special.expn(1, tau[i]-tau[0:i]) , tau[0:i]) + \
               2*np.pi*np.trapz(S[i+1:]*scipy.special.expn(1, tau[i+1:]-tau[i]), tau[i+1:])

        #######################################################################
    
        # Evaluate q
        q[i] = 2*J[0]*scipy.special.expn(3,tau[i]) - \
               2*J[1]*scipy.special.expn(3,tau[-1]-tau[i] ) + \
               2*np.pi*np.trapz(S[0:i+1]*scipy.special.expn(2, tau[i]-tau[0:i+1]) , tau[0:i+1]) - \
               2*np.pi*np.trapz(S[i:]*scipy.special.expn(2, tau[i:]-tau[i])       , tau[i:]) 

    # Update the source term
    if scattering:
        S = np.multiply((1-omega),I_b) + np.multiply(omega/(4*np.pi),G)
    else:
        S = I_b
   
    # Find the L2 norm of the error the old and the new geuess for G
    res = np.linalg.norm(G - G_old)
    print('res = {}'.format(res))
    
    # Make the current guess the old guess
    G_old = copy.deepcopy(G)

# save the x, G, and q arrays to a compressed numpy data file ( .npz file)
np.savez(output_file_name, x_analytical=x, G_analytical=G, q_analytical=q)

plt.figure()
plt.plot(x, T, '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$T \ \left[ K \right]$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure()
plt.plot(x, sigma_a, '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$\sigma_a \ \left[ \frac{1}{m} \right]$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure()
plt.plot(x, sigma_s, '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$\sigma_s \ \left[ \frac{1}{m} \right]$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure()
plt.plot(x, sigma_e, '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$\sigma_e \ \left[ \frac{1}{m} \right]$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure()
plt.plot(x, omega, '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$\omega$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure()
plt.plot(x, tau, '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$\tau$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure()
plt.plot(x, G , '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$G \ \left[ \frac{W}{m^2 \cdot Sr} \right]$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure()
plt.plot(x, q , '-k', label=r'Analytical')
plt.xlabel(r'$x$', fontsize = 20)
plt.ylabel(r'$q \ \left[ \frac{W}{m^2} \right]$', fontsize = 20)
plt.legend(loc = 'best')

plt.show()
