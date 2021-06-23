import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import test_gas

"""
This is the PFR veriosn of "Auto_reactor_modell_surf.py". 
The difference between this two is a way of simulation. In this PFR, simulation is conducted on 
every single reactor individually till the steady state of every reactor is reached.
Next, the new iteration starts which brings at the end the network of reactors to steady state.
This is in contrary to the previous version where simulation is always conducted on all reactors simultaneously. 

- working & by iterations ~ 100 delivers good results
- much faster than previous version
- it is possible to simulate any number of reactors (much faster)

TO DO:
- Implement mech_13.yaml or other mechanisms as gri30.yaml
    It has to have a 'Pt_surf' mechanism in the reaction to simulate it. There can be no reaction as 
    surface reaction simulated that has no surface reaction mechanism in it. 

- Implement the reduced mechanism
    Under which conditions the machanism should be reduced? The 'methane_pox_on_pt.cti' is not working because 
    only the gas part of the reaction is taken into simulation. It has to be corrected!
    Temp, Press, X, time, k?

- Carbon (H2) deposition mechanism

NOT WORKING:
- Removal of one component from the gas mixture after every reactor to simulate deposition
"""
#######################################################################
# Input Parameters
#######################################################################

# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.
# IT MUST DIVISBLE BY 3!!!!
n_steps = int(input('Number of reactors in the pfr system. The number must be divisible by 4!\n'))
iters = int(input('Number of iterations for the pfr-reactor cylcle?\n'))
mechanism = int(input('Press 0 for automatic gri_30, 1 to choose reduced gri30 or Press 2 to choose mech_13\n'))

# starts counting time
start_time = time.time()

# unit conversion factors to SI
cm = 0.01
minute = 60.0

# inlet temperature [K]
T_0 = 373.0
# wall temperature [K]
T_wall = 873.0
# constant pressure [Pa]
pressure = ct.one_atm
composition_0 = 'C3H8:10, H2:1'
initial_state = T_0, pressure, composition_0

# Reactor geometrie
vol_flow_rate = 1.6664e-6 # [m3/s]
length = 0.1  # *approximate* PFR length [m]
area = 0.00024  # cross-sectional area [m**2]
u_0 = vol_flow_rate/area # inflow velocity [m/s]

if mechanism > 0:
    # Calling a function from another script to choose the desired mechanism after reduction
    gas, reaction_mechanism1 = test_gas.gas_funct(mechanism, initial_state)
else:
    reaction_mechanism1 = 'gri30.yaml'
    # import the gas model and set the initial conditions
    gas = ct.Solution(reaction_mechanism1)

gas.TPX = initial_state

mass_flow_rate = u_0 * gas.density * area
dz = length / n_steps
r_vol = area * dz


# Create ideal gas reactors.
reactors = [ct.IdealGasReactor(contents = gas, energy = 'on') for i in range(n_steps)]

# Add a volume to reactors
for i, reactors_obj in enumerate(reactors):
    reactors_obj.volume = r_vol

# create reservoirs for heat exchange with gas2
gas2 = ct.Solution('liquidvapor.yaml', 'water')
gas2.TPX = 873.0, 5e+5, 'H2O:1'
gas2()
n_walls = n_steps // 2
heat_reserv = [ct.Reservoir(gas2) for i in range(n_walls)]

# # set temperature in these Reservoirs
# for reserv in heat_reserv:
#     reserv.thermo.TPX = T_wall, pressure, composition_0

# create wall for heat exchange. What Wall area [m^2] and Overall heat transfer coefficient [W/m^2]??
for i in range(n_walls):
    ct.Wall(heat_reserv[i], reactors[i + n_steps // 4], A=dz*0.006, U=1.2)

# create a reservoir to represent the reactor immediately upstream
upstreams = [ct.Reservoir(gas) for i in range(n_steps + 1)]

# create a downstream reservoir
downstreams = [ct.Reservoir(gas) for i in range(n_steps)]


# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
mflows = []
for i in range(n_steps):
    mflows.append(ct.MassFlowController(upstreams[i], reactors[i], mdot=mass_flow_rate))

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference. What K?!
vpress = []
for i in range(n_steps):
    vpress.append(ct.PressureController(reactors[i], downstreams[i], master=mflows[i], K=1e-5))
    # vpress.append(ct.Valve(reactors[i], downstreams[i], K=1e-5))


# "Länge" des Gesamtreaktors nach jedem Teilreaktor" ??
z_vec = []
for i in range(n_steps):
    z_vec.append(i * dz)

# Zustände in den Reaktoren
dict_states = {}
for i in range(n_steps):
    dict_states[i] = ct.SolutionArray(reactors[i].thermo)

# Create an empty dictionary to store Temperatures of individual reactors
dict_T = {}
# Adding arrays to this dictionary
for i in range(n_steps):
    dict_T[i] = []

# Create Lists for evaluation
temp_profile = []
p_profile = []
H2_profile = []
X_CH4_profile = []
X_O2_profile = []
X_CO2_profile = []
X_H2O_profile = []

# Create a list for velocity for every reactor
u = []
# For residence time in each reactor
t_r = []
# Whole time
t = []
# Density in each reactor
r_dens = []

q1 = ct.Quantity(gas)
print("There is ", q1.moles, "moles gas in system\n")

# The outer loop iterate through the given number of iterations for every pfr cycle
for n in range(iters):
    # This loop iterate through all reactors
    for i in range(n_steps):

        # print(i, " Reactor")
        # Define sim as simulation in one reactor
        sim = ct.ReactorNet([reactors[i]])

        # Set the state of the reservoir to match that of the previous reactor
        if i > 0:
            gas.TDY = upstreams[i].thermo.TDY
        else:
            gas.TPX = T_0, pressure, composition_0

        # integrate the i reactor forward in time until steady state is reached

        sim.reinitialize()
        sim.advance_to_steady_state()

        downstreams[i].syncState()
        # print("Upstream: ", i, upstreams[i].thermo.X[0], "\n")
        # print("Reactor: ", i, reactors[i].thermo.X[0], "\n")
        # print("Downstream: ", i, downstreams[i].thermo.X[0], "\n")

        # # Save temperature & thermo.state in each Reaktor
        # for j in range(n_steps):
        #     dict_T[j].append(reactors[j].T)
        #     dict_states[j].append(reactors[j].thermo.state)

        upstreams[i + 1].thermo.TPX = downstreams[i].thermo.TPX
        upstreams[i + 1].syncState()

    print('Reaktordurchlauf:\t', n + 1)

# Creating lists with temp, pressure and species conc. for diagrams.
for k in range(n_steps):
    # compute velocity and transform into time
    u.append(mass_flow_rate / area / reactors[k].density)
    t_r.append(reactors[k].mass / mass_flow_rate)  # residence time in this reactor
    t.append(t_r[k])
    r_dens.append(reactors[k].thermo.density)

    temp_profile.append(reactors[k].T)
    p_profile.append(reactors[k].thermo.P)
    H2_profile.append(reactors[k].thermo.X[gas.species_index('H2')])
    X_CH4_profile.append(reactors[k].thermo.X[gas.species_index('C3H8')])
    X_O2_profile.append(reactors[k].thermo.X[gas.species_index('O2')])
    X_CO2_profile.append(reactors[k].thermo.X[gas.species_index('CO2')])
    X_H2O_profile.append(reactors[k].thermo.X[gas.species_index('H2O')])

# At the outlet from Reactor
gas.TDY = reactors[n_steps - 1].thermo.TDY
gas()

# Plots:
# plt.figure()
# plt.subplot(231)
# plt.plot(z_vec, temp_profile)
# plt.subplot(232)
# plt.plot(z_vec, p_profile)
# plt.subplot(233)
# plt.plot(z_vec, X_CH4_profile)
# plt.subplot(234)
# plt.plot(z_vec, X_O2_profile)
# plt.subplot(235)
# plt.plot(z_vec, X_CO2_profile)
# plt.subplot(236)
# plt.plot(z_vec, X_H2O_profile)
#
# plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Plots for evaluation of reactor')
ax1.plot(z_vec, temp_profile)
ax1.set_ylabel('Temperature profile')

ax2.plot(z_vec, X_CH4_profile)
ax2.set_ylabel('C3H8 profile')

ax3.plot(z_vec, H2_profile)
ax3.set_ylabel('H2 profile')
ax2.set_xlabel('Reactor length')

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))