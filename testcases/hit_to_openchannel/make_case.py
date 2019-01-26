#!/usr/bin/env python3

import argparse
import csv
import json
import numpy as np
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--debug', action='store_true')
parser.add_argument('base_json', type=argparse.FileType('r'))
parser.add_argument('Re_W', type=float)
parser.add_argument('St_eta', type=float)
parser.add_argument('MLR', type=float)
parser.add_argument('R', type=float)
parser.add_argument('chi', type=float)
parser.add_argument('C', type=float)
args = parser.parse_args()

mc = json.load(args.base_json)

def output(*strings):
    if args.debug:
        print(*strings)

#### Fixed values #################################################################################

## Inflow-outflow length
L = 0.16				# [m]

## Inflow-outflow width
W = 0.04				# [m]

## Inlet temperature
T_0 = 300.0				# [K]

## Maximum expected temperature
T_max = 1000.0				# [K]

## Inlet fluid density
rho_0 = 1.2				# [kg/m^3]

## Fluid dynamic viscosity (constant)
mu = 1.9e-5				# [Pa s]

## Fluid isobaric heat capacity
C_p = 1012.0				# [J/(kg K)]

## Fluid isochoric heat capacity
C_v = 723.0				# [J/(kg K)]

## Particle Nusselt number
Nu = 2.0				# [-]

## Particle density
rho_p = 8900.0				# [kg/m^3]

## Particle emissivity (absortivity)
epsilon_p = 0.4				# [-]

## Particle absorption and scattering efficiencies
Q_a = 0.4				# [-]
Q_s = 0.7				# [-]

## Mean rms velocity fluctuation (HIT) to inlet velocity ratio
u_rms_U_0 = 0.12			# [-]

## Re_W to Re_lambda conversion coefficient for this flow
Re_W_Re_lambda = ( 43.0**2.0 )/5217.0	# [-]

## Number of grid points
N_x = mc['configs'][1]['Grid']['xNum']
N_y = mc['configs'][1]['Grid']['yNum']
N_z = mc['configs'][1]['Grid']['zNum']

## CFL limit
CFL = 0.25

## Factor of Safety ratio of scales
FoS = 1.0/5.0

#### Input values #################################################################################

## Bulk Reynolds number
Re_W = args.Re_W

## Kolmogorov Stokes number
St_eta = args.St_eta

## Mass loading ratio
MLR = args.MLR

## Dimensionless radiation intensity
R = args.R

## Heat capacity ratio
chi = args.chi

## Radiation to convection heat transfer ratio
C = args.C

#### Output values ################################################################################

output()

## Inlet velocity
U_0 = ( mu*Re_W )/( W*rho_0 )								# [m/s]
output( 'Inlet velocity: \t \t  U_0       =', U_0, '\t     [m/s]' )

## Taylor Reynolds number (HIT)
Re_lambda = np.sqrt( Re_W_Re_lambda*Re_W )						# [-]
#output( 'Taylor Reynolds number: \t  Re_lambda =', Re_lambda, '\t \t     [-]' )

## Mean rms velocity fluctuation (HIT)
u_rms = U_0*u_rms_U_0									# [m/s]
#output( 'Mean velocity fluctuation: \t  u_rms     =', u_rms, '   [m/s]' )

## Taylor length scale (HIT)
lambda_ = ( Re_lambda*mu )/( rho_0*u_rms )						# [m]
#output( 'Taylor length scale: \t          lambda    =', lambda_, '      [m]' )

## Mean KTE dissipation rate (HIT)
epsilon = ( 15.0*mu*( u_rms**2.0 ) )/( rho_0*( lambda_**2.0 ) )				# [m^2/s^3]
#output( 'Mean TKE dissipation rate: \t  epsilon   =', epsilon, '          [m^2/s^3]' )

## Kolmogorov length scale (HIT)
eta = ( ( mu**3.0 )/( ( rho_0**3.0 )*epsilon ) )**( 1.0/4.0 )				# [m]
#output( 'Kolmogorov length scale: \t  eta       =', eta, '     [m]' )

## Kolmogorov time scale
tau_eta = np.sqrt( mu/( rho_0*epsilon ) )						# [s]
#output( 'Kolmogorov time scale: \t \t  tau_eta   =', tau_eta, '\t     [s]' )

## Particle diameter
d_p = np.sqrt( ( 18.0*mu*St_eta*tau_eta )/( rho_p ) )					# [m]
output( 'Particle diameter: \t \t  d_p       =', d_p, '     [m]' )

## Inflow-outflow domain volume
V = L*W*W										# [m^3]
#output( 'Inflow-outflow domain volume: \t  V         =', V, '[m^3]' )

## Particle mass
m_p = ( np.pi/6.0 )*rho_p*( d_p**3.0 )							# [kg]
#output( 'Particle mass: \t \t \t  m_p       =', m_p, '     [kg]' )

## Inflow-outflow number of particles
N_p = int( ( rho_0*V*MLR )/m_p )							# [-]
output( 'I-O number of particles: \t  N_p       =', N_p, '\t     [-]' )

## Particle number density
n_0 = N_p/V										# [-]
#output( 'Particle number density: \t  n_0       =', n_0, '         [-]' )

## HIT number of particles
N_p_HIT = int( n_0*( W*W*W ) )								# [-]
output( 'HIT number of particles: \t  N_p       =', N_p_HIT, '\t     [-]' )

## Radiation intensity
I_0 = ( 4.0*T_0*rho_0*C_p*U_0*C )/( np.pi*epsilon_p*n_0*L*( d_p**2.0 ) )		# [W/m^2]
output( 'Total Radiation intensity: \t  I_0       =', I_0, '\t     [W/m^2]' )
output( '1/2 Radiation intensity: \t  1/2 I_0   =', 0.5*I_0, '\t     [W/m^2]' )

## Particle isochoric heat capacity
C_v_p = ( C_p*chi )/MLR									# [J/(kg K)]
output( 'Particle isochoric heat capacity: C_v_p     =', C_v_p, '     [J/(kg K)]' )

## Fluid thermal conductivity
kappa = ( epsilon_p*I_0*d_p )/( 4.0*Nu*T_0*R )						# [W/(m K)]
#output( 'Thermal conductivity: \t \t  kappa     =', kappa, '\t     [W/(m K)]' )

## Prandtl number
Pr = ( C_p*mu )/kappa									# [-]
output( 'Prandtl number: \t \t  Pr        =', Pr, '\t     [-]' )

## Particle convection coefficient
h = ( Nu*kappa )/d_p									# [W/(m^2 K)]
output( 'Particle convection coefficient:  h         =', h, '\t     [W/(m^2 K)]')

## Large-eddy length scale (HIT)
l_0 = 0.19*W										# [m]
#output( 'Large-eddy length scale (HIT): \t  l_0 \t    =', l_0, '\t \t     [m]' )

## Large-eddy velocity (HIT)
u_0 = ( ( mu/rho_0 )*( Re_lambda**2.0 ) )/( 15.0*l_0 )					# [m/s]
#output( 'Large-eddy velocity (HIT): \t  u_0 \t    =', u_0, '\t     [m/s]' )

## Large-eddy time scale (HIT)
t_0 = l_0/u_0										# [s]
output( 'Large-eddy time scale (HIT): \t  t_0 \t    =', t_0, '\t     [s]' )

## Turbulent kinetic energy (HIT)
k_0 = ( 3.0/2.0 )*( u_0**2.0 )								# [m^2/s^2]
output( 'Turbulent kinetic energy (HIT):   k_0 \t    =', k_0, '       [m^2/s^2]' )

## Speed of sound (air)
#c = 331.5 + 0.6*( T_0 - 273.15 )							# [m/s]
c = 331.5 + 0.6*( T_max - 273.15 )							# [m/s]
#output( 'Speed of sound (air):  \t \t  c \t    =', c, ' \t \t     [m/s]' )

## Flow acoustic time step
delta_t_c = CFL*( L/N_x )/c								# [s]
output( 'Flow acoustic time step:  \t  delta_t_c =', delta_t_c, ' [s]' )

## Particle relaxation time: time required for a particle to adapt to the flow eddies
tau_p = St_eta*tau_eta									# [s]
#output( 'Particle relaxation time: \t  tau_p     =', tau_p, '\t     [s]' )

## Particle thermal relaxation time: time required for the thermal energy of a particle to adapt to the flow
tau_th = ( 3.0*tau_p*C_v_p*Pr )/( Nu*C_p )						# [s]
#output( 'Particle thermal relaxation time: tau_th    =', tau_th, '\t     [s]' )

## Optical thickness
tau = ( ( np.pi*d_p*d_p )/4.0 )*W*n_0							# [-]
#output( 'Optical thickness: \t \t  tau \t    =', tau, '\t     [-]' )

## Shadow fraction in the limit of small particles in comparison to the duct width
SF = 1.0 - np.exp( ( -1.0 )*tau )							# [-]
#output( 'Shadow fraction: \t \t  SF \t    =', SF, '\t     [-]' )

## Radiation albedo
omega = Q_s/( Q_a + Q_s )								# [-]
#output( 'Radiation albedo: \t \t  omega     =', omega, '    [-]' )

## Radiation time scale
tau_rad = ( rho_0*C_p*T_0*W )/( I_0*SF*( 1 - omega ) )					# [s]
#output( 'Radiation time scale: \t \t  tau_rad   =', tau_rad, '\t     [s]' )

## Ratio between minimal tau scales and acoustic time step
ratio_delta_t = int( FoS*( min( tau_eta, tau_p, tau_th, tau_rad )/delta_t_c ) )		# [-]
output( 'Ratio scales to c time step: \t  ratio     =', ratio_delta_t, '\t \t     [-]' )

## Frequency copy between domains in number of flow time step
ratio_copy = int( ( ( L/N_x )/U_0 )/delta_t_c )						# [-]
output( 'Frequency copy between domains:   # delta_t =', ratio_copy, '\t \t     [-]' )

output()

#### Fill in json config ##########################################################################

if not args.debug:
    # Round up number of particles to fit tiling
    tiles_0 = (mc['configs'][0]['Mapping']['tiles'][0] *
               mc['configs'][0]['Mapping']['tiles'][1] *
               mc['configs'][0]['Mapping']['tiles'][2])
    if N_p_HIT % tiles_0 > 0:
        N_p_HIT += tiles_0 - (N_p_HIT % tiles_0)
    tiles_1 = (mc['configs'][1]['Mapping']['tiles'][0] *
               mc['configs'][1]['Mapping']['tiles'][1] *
               mc['configs'][1]['Mapping']['tiles'][2])
    if N_p % tiles_1 > 0:
        N_p += tiles_1 - (N_p % tiles_1)
    # Adjust stagger factor to be compatible with copy frequency
    ratio_delta_t = min(ratio_delta_t, ratio_copy)
    factor = 1
    while ratio_copy % factor != 0 or ratio_copy // factor > ratio_delta_t:
        factor += 1
    ratio_delta_t = ratio_copy // factor
    # Fill in variable parameters
    mc['configs'][0]['Integrator']['fixedDeltaTime'] = delta_t_c
    mc['configs'][0]['Flow']['prandtl'] = Pr
    mc['configs'][0]['Flow']['turbForcing']['t_o'] = t_0
    mc['configs'][0]['Flow']['turbForcing']['K_o'] = k_0
    mc['configs'][0]['Particles']['initNum'] = N_p_HIT
    mc['configs'][0]['Particles']['maxNum'] = N_p_HIT
    mc['configs'][0]['Particles']['convectiveCoeff'] = h
    mc['configs'][0]['Particles']['heatCapacity'] = C_v_p
    mc['configs'][0]['Particles']['diameterMean'] = d_p
    mc['configs'][0]['Particles']['staggerFactor'] = ratio_delta_t
    mc['configs'][1]['BC']['xBCLeftInflowProfile']['addedVelocity'] = U_0
    mc['configs'][1]['Integrator']['fixedDeltaTime'] = delta_t_c
    mc['configs'][1]['Flow']['prandtl'] = Pr
    mc['configs'][1]['Flow']['initParams'][2] = U_0
    mc['configs'][1]['Particles']['maxNum'] = int(N_p * 1.2)
    mc['configs'][1]['Particles']['convectiveCoeff'] = h
    mc['configs'][1]['Particles']['heatCapacity'] = C_v_p
    mc['configs'][1]['Particles']['diameterMean'] = d_p
    mc['configs'][1]['Particles']['feeding']['addedVelocity'][0] = U_0
    mc['configs'][1]['Particles']['staggerFactor'] = ratio_delta_t
    if mc['configs'][1]['Radiation']['type'] == 'DOM':
        mc['configs'][1]['Radiation']['yHiIntensity'] = I_0/2
        mc['configs'][1]['Radiation']['yLoIntensity'] = I_0/2
    elif mc['configs'][1]['Radiation']['type'] == 'Algebraic':
        mc['configs'][1]['Radiation']['absorptivity'] = Q_a
        mc['configs'][1]['Radiation']['intensity'] = I_0
    else:
        assert false
    mc['copyEveryTimeSteps'] = ratio_copy
    # Dump final json config
    json.dump(mc, sys.stdout, indent=4)
