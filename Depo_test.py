import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import V_N
import Reduction_mech
import json

# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.

#######################################################################
# Input Parameters
#######################################################################
n_steps = int(input('Number of reactors in the pfr system. The number must be divisible by 4!\n'))
# iters = int(input('Number of iterations for the pfr-reactor cylcle?\n'))
mechanism = int(input('Press 0 for automatic gri_30, 1 to choose reduced gri30 or Press 2 to choose mech_13\n'))
remove = int(input('Should all heavy components be removed after every reactor? 1- yes 0- no\n'))

# starts counting time
start_time = time.time()

# DATA
T_0_n = 273.0
# inlet temperature [K]
T_0 = 473.0
# wall temperature [K]
T_wall = 873.0
# constant pressure [Pa]
pressure = ct.one_atm
composition_0 = 'C3H8:10, H2:1'
initial_state = T_0, pressure, composition_0
reactive_state = T_wall, pressure, composition_0

# Here the reaction mechanism will be choosed
if mechanism > 0:
    # Calling a function from another script to choose the desired mechanism after reduction
    gas, reaction_mechanism1 = Reduction_mech.gas_funct(mechanism, reactive_state)
else:
    reaction_mechanism1 = 'gri30.yaml'
    # import the gas model and set the initial conditions
    gas = ct.Solution(reaction_mechanism1)

# State with standard temp 273 for vol. calculation
gas.TPX = T_0_n, pressure, composition_0
gas()

# Reactor geometry
vol_flow_rate_N = 1.6664e-6  # [m3/s] volumen flow standarized for 273K
vol_flow_rate = V_N.vol_norm(vol_flow_rate_N, gas, T_0, pressure)  # [m3/s]
length = 0.1  # *approximate* PFR length [m]
area = 0.00024  # cross-sectional area [m**2]
u_0 = vol_flow_rate / area  # inflow velocity [m/s]
height = 0.006  # [m]
depth = 0.04  # [m]

# Back to initial state
gas.TPX = initial_state
gas()

# Calculate other data
mass_flow_rate_0 = u_0 * gas.density * area
dz = length / n_steps
r_vol = area * dz
A_wall = 2 * height * dz + 2 * depth * dz

# create a new reactor
reactor = ct.IdealGasReactor(contents=gas, energy='on')
reactor.volume = r_vol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream = ct.Reservoir(gas, name='upstream')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream = ct.Reservoir(gas, name='downstream')

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
m = ct.MassFlowController(upstream, reactor)
m.mass_flow_rate = mass_flow_rate_0

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v = ct.PressureController(reactor, downstream, master=m, K=1e-6)

# create reservoirs for heat exchange with gas2 and for enviroment temp gas 3
gas_600 = ct.Solution('liquidvapor.yaml', 'water')
gas_600.TPX = 873.0, 5e+5, 'H2O:1'
# gas2()
heat_reserv = ct.Reservoir(gas_600)

gas_200 = gas_600
gas_200.TPX = 473.0, 5e+5, 'H2O:1'
# gas3()
env_reserv = ct.Reservoir(gas_200)

####### Initializing some lists for holding results ##########
# "Länge" des Gesamtreaktors nach jedem Teilreaktor" ??
z_vec = []
for i in range(n_steps + 1):
    z_vec.append(i * dz)
# moles
moles = [0] * (n_steps + 1)
# mass flows
mflow = [0] * (n_steps + 1)
mflow[0] = mass_flow_rate_0
# mole flow
mol_flow = [0] * (n_steps + 1)
mol_flow[0] = mass_flow_rate_0 / reactor.thermo.mean_molecular_weight
# velocity through the reactor
u = [0] * (n_steps + 1)
u[0] = u_0
# This dictionary is for saving the next reactor state. After every loop it will be reconstructed from this data.
state_list = [0] * (n_steps + 1)
state_list[0] = reactor.thermo.TPY
# mass in every reactor
mass = [0] * (n_steps + 1)
mass[0] = reactor.mass
# pressure in every reactor
p_profile = [0] * (n_steps + 1)
p_profile[0] = reactor.thermo.P
# temp profile
T_profile = [0] * (n_steps + 1)
T_profile[0] = reactor.T
# density profile
D_profile = [0] * (n_steps + 1)
D_profile[0] = reactor.density
# Time in reactive part
t_react = 0

# Create a dictionary for storing deposition components and their moles
depo = {}
depo_r = {}
mass_depo = [0] * (n_steps + 1)

# iterate through the PFR cells

"""
An idea on how to solve the variable flow rate problem:
1) Create data structures to store each individual reactor state. 
2) After each reactor (n_step) store the new reactor state and the new mass flow at the end.
3) The gas thermodynamics and mass flow of each subsequent reactor 
    will be defined by the previous state stored in the corresponding structure.
4) Simulation.
5) Removal of the compound. Calculation of a new mass flow based on the previous one.
6) The new iteration may overwrite the previous one. (Maybe??). 
7) Create mass flow graphs through the reactor network.

Problems: 
8) Link mass flow to reactor pressure!!! It seems that the pressure is held constant in the reactor. 
 -> Pressure should remain const.
9) What should be constant, mass flow or e.g. volume flow, i.e. velocity only?!(Key!!!)
 -> mass flow rate is also const.
10) Compare the two simulation methods without deposition. 
11) Write the deposition in such a way that successive elements are removed "additively".-> Done
12) How to calculate density with deposition? !!!!!!!! -> deprecated
13) Create a dict that records and ideally adds deposits from each reactor -> Done
14) Check which heaviest elements can be deposited in gri30 to make it faster than in mech_13 -> Done (Only ~ 40 g/mol)
15) Plot deposition through the length of reactor 
"""

# for i in range (iters) :
#     # To reset the old state
#     reactor = ct.IdealGasReactor(contents = gas, energy = 'on')
for n in range(n_steps):
    # delete all existing walls and previos state in the reactor
    # reactor.thermo.TDY.clear()
    reactor.walls.clear()
    # create wall for heat exchange. What Wall area [m^2] and Overall heat transfer coefficient [W/m^2]??
    if n < (3 * (n_steps // 4)):
        ct.Wall(heat_reserv, reactor, A=A_wall, U=50)
        t_react += reactor.mass / mflow[n]
    else:
        ct.Wall(env_reserv, reactor, A=A_wall, U=50)

    sim2 = ct.ReactorNet([reactor])

    # Set the state of the reservoir to match that of the previous reactor
    gas.TPY = state_list[n]
    upstream.syncState()
    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()
    sim2.advance_to_steady_state()

    T = reactor.T
    P = reactor.thermo.P
    Y = reactor.Y

    p_profile[n + 1] = reactor.thermo.P
    T_profile[n + 1] = T

    ############## REMOVE ONE COMPONENT FROM THE REACTOR NET ###########

    # remove water after every reactor
    if remove == 1:

        # These are entry data if no compounds was removed
        mass[n + 1] = reactor.mass
        mass_flow_now = mflow[n]
        mflow[n + 1] = mass_flow_now
        mol_flow[n + 1] = mass_flow_now / reactor.thermo.mean_molecular_weight
        state_list[n + 1] = (T, P, Y)
        mass_depo[n] = 0.0

        # For every compound in the reactor show only these with mol. weight ...
        # These will be set to 0.0
        for j, compound in enumerate(reactor.thermo.species()):

            if 30.0 > reactor.thermo.molecular_weights[j] > 20.00:

                name_Y_excl = compound.name
                Y_excl = gas.species_index(compound.name)

                if name_Y_excl not in depo.keys():
                    depo[name_Y_excl] = 0.0
                    depo_r[name_Y_excl] = []

                # Mass in reactor
                mass[n + 1] -= reactor.mass * Y[Y_excl]

                # Density in reactor after deposition
                # D -= reactor.mass * Y[Y_excl] / reactor.volume # Is that correct?!

                # Part of the mass flow that is deposited in one slicereactor
                mass_flow_rate_down = mass_flow_now * Y[Y_excl]
                # Sum of deposited mass
                mass_depo[n] += mass_flow_rate_down
                # Set new mass flow rate in upstream
                mflow[n + 1] -= mass_flow_rate_down
                # Set new mol flow rate in the upstream
                mol_flow[n + 1] -= mass_flow_rate_down / gas.molecular_weights[Y_excl]

                # Set new mass fraction after deposition (later it can be a function)
                Y[Y_excl] = 0.0
                state_list[n + 1] = (T, P, Y)

                # Save deposited compound and moles in a dict.
                depo[name_Y_excl] += mass_flow_rate_down
                depo_r[name_Y_excl].append(mass_flow_rate_down)

        # Set a new value for mass flow in the upstream
        m.mass_flow_rate = mflow[n + 1]
        u[n+1] = mflow[n+1]/area*reactor.density
        upstream.syncState()
    else:
        mass[n + 1] = reactor.mass
        mflow[n + 1] = mflow[n]
        u[n + 1] = mflow[n + 1] / area * reactor.density
        mol_flow[n + 1] = mflow[n + 1] / reactor.thermo.mean_molecular_weight
        state_list[n + 1] = (T, P, Y)
        upstream.syncState()
########################################################################
gas()

# Creating a file to store deposited compounds and their mass flow rate
if remove == 1:
    f = open("Deposition.txt", "w")
    f.write("Comp.: MFR[kg/s]:\n\n")
    for key in depo:
        f.write(key + ": " + str(depo[key]) + "\n")
    f.close()

# Look for the compound that has a biggest sum of deposition
# Compound with biggest deposition
depo_comp_name = min(depo,key=depo.get)
print(depo_comp_name)
depo_comp_reactor = []
depo_comp_reactor.append(0.0)
for i in range(n_steps):
    depo_comp_reactor.append(depo_r[depo_comp_name][i])

###################################################################
# Compare Results in matplotlib
###################################################################
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Plots for evaluation of reactor T & p')
ax1.plot(z_vec, p_profile)  # [s[1] for s in state_list]
ax1.set_ylabel('Pressure [Pa]')

ax2.plot(z_vec, T_profile)
ax2.set_ylabel('Temperature [K]')

ax3.plot(z_vec, u)
ax3.set_ylabel('Velocity [m/s]')

ax2.set_xlabel('Reactor length [m]')

fig2, (ax4, ax5, ax6) = plt.subplots(3, 1)
fig2.suptitle('Flow rates & mass across the reactor')
ax4.plot(z_vec, mflow)
ax4.set_ylabel('mass flow rate [kg/s]')

ax5.plot(z_vec, mol_flow)
ax5.set_ylabel('mol flow rate [mol/s]')

ax6.plot(z_vec, mass)
ax6.set_ylabel('mass across reactor [kg]')

ax6.set_xlabel('Reactor length [m]')

fig3, (ax7, ax8) = plt.subplots(2,1)
fig3.suptitle('Deposition across the reactor')

ax7.plot(z_vec, mass_depo)
ax7.set_ylabel('sum of depo across reactor [kg/s]')

ax8.plot(z_vec, depo_comp_reactor)
ax8.set_ylabel(depo_comp_name + ' depo across the reactor')

ax8.set_xlabel('Reactor length [m]')

plt.show()

print("\nResidence time in reactive Part = ", t_react, " s")

print("\n--- %s seconds ---" % (time.time() - start_time))
