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
n_steps = int(input('Number of reactors in the pfr system. The number must be divisible by 3!\n'))
iters = int(input('Number of iterations for the pfr-reactor cylcle?\n'))
walls = int(input('Press 1 to activate walls for heat transfer\n'))
mechanism = int(input('Press 0 for automatic methane, 1 to choose reduced methane pox on pt or Press 2 to choose mech_13\n'))

# starts counting time
start_time = time.time()

# unit conversion factors to SI
cm = 0.01
minute = 60.0

T_0 = 1000.0  # inlet temperature [K]
pressure = ct.one_atm  # constant pressure [Pa]
composition_0 = 'CH4:1, O2:1.5, AR:0.1'
# oder 'C3H8:10, H2:1'
initial_state = T_0, pressure, composition_0

# catalyst
# length = 0.3 * cm  # Catalyst bed length
# area = 1.0 * cm**2  # Catalyst bed area
cat_area_per_vol = 1000.0 / cm  # Catalyst particle surface area per unit volume
initial_coverage_catalyst = 'O(S):0.00, PT(S):0.01, H(S):0.99'  # Zeile selbst hizugefügt. Wert aus Mechanismus

# Warum solche Angaben?????????????
length = 1.5e-7  # *approximate* PFR length [m]
u_0 = .006  # inflow velocity [m/s]
area = 1.e-4  # cross-sectional area [m**2]

if mechanism > 0:
    # Calling a function from another script to choose the desired mechanism after reduction
    gas, reaction_mechanism1 = test_gas.gas_funct(mechanism, initial_state)
else:
    reaction_mechanism1 = 'methane_pox_on_pt.cti'
    # import the gas model and set the initial conditions
    gas = ct.Solution(reaction_mechanism1, 'gas')


gas.TPX = initial_state

mass_flow_rate = u_0 * gas.density * area
dz = length / n_steps
r_vol = area * dz

# Automated way of creating surfaces is applied
surfs = [ct.Interface(reaction_mechanism1, 'Pt_surf', [gas]) for i in range(n_steps)]

for i, surf_obj in enumerate(surfs):
    surf_obj.TP = T_0, ct.one_atm
    # print(surf_obj.TP)

# catalyst area in one reactor
cat_area = cat_area_per_vol * r_vol

# Create ideal gas reactors
reactors = [ct.IdealGasReactor(gas) for i in range(n_steps)]

# Add a volume to reactors
for i, reactors_obj in enumerate(reactors):
    reactors_obj.volume = r_vol

# Add walls between reactors on the opposite side of system
if walls == 1:
    n_walls = n_steps // 3
    # print(n_walls)

    for i in range(n_walls):
        ct.Wall(reactors[i], reactors[n_steps - 1 - i], A=0.001, U=1.2)
    print('Walls are active!')
else:
    print('Wall are not active!')

# create a reservoir to represent the reactor immediately upstream
upstreams = [ct.Reservoir(gas) for i in range(n_steps)]

# create a downstream reservoir
downstreams = [ct.Reservoir(gas) for i in range(n_steps)]

# Add the reacting surface to the reactor. The area is set to the desired
# catalyst area in the reactor.
rsurfs = []
# Between 1/3 n_steps & 2/3 n_steps
for i in range(n_steps // 3, 2 * (n_steps // 3)):
    rsurfs.append(ct.ReactorSurface(surfs[i], reactors[i], A=cat_area))

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
mflows = []
for i in range(n_steps):
    mflows.append(ct.MassFlowController(upstreams[i], reactors[i], mdot=mass_flow_rate))

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
vpress = []
for i in range(n_steps):
    vpress.append(ct.PressureController(reactors[i], downstreams[i], master=mflows[i], K=1e-5))


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

# The outer loop iterate through the given number of iterations for every pfr cycle
for n in range(iters):
    # This loop iterate through all reactors
    for i in range(n_steps):

        # print(i, " Reactor")
        # Define sim as simulation in one reactor
        sim = ct.ReactorNet([reactors[i]])

        # Set the state of the reservoir to match that of the previous reactor
        if i > 0:
            gas.TDY = reactors[i - 1].thermo.TDY
        else:
            gas.TPX = T_0, pressure, composition_0

        upstreams[i].syncState()

        # integrate the i reactor forward in time until steady state is reached

        sim.reinitialize()
        sim.advance_to_steady_state()

        # print("Upstream: ", i, upstreams[i].thermo.X[0], "\n")
        # print("Reactor: ", i, reactors[i].thermo.X[0], "\n")
        # print("Downstream: ", i, downstreams[i].thermo.X[0], "\n")

        # # Save temperature & thermo.state in each Reaktor
        # for j in range(n_steps):
        #     dict_T[j].append(reactors[j].T)
        #     dict_states[j].append(reactors[j].thermo.state)

    ############# This is an attempt to remove one component from the solution #############
    #
    # #     What species are present in the system
    # #     print(gas.species_names)
    # #     print(gas.X)
    # #     for ii, i_r in enumerate(reactors):
    # #         print(ii, i_r.thermo.X[0])
    #     # This is not correct, it has to be done by relative reset H2.
    #
    #     # print(reactors[i].thermo.molecular_weights)
    #     for j, compound in enumerate(reactors[i].thermo.species()):
    #         if reactors[i].thermo.molecular_weights[j] < 4.0:
    #             # print(compound)
    #             X_excl = gas.species_index(compound.name)
    #             T,P,X = reactors[i].thermo.TPX
    #             X[X_excl] = X[X_excl]*0.1
    #             reactors[i].thermo.TPX = T,P,X
    #             reactors[i].syncState()
    #             # print(reactors[i].thermo.X[X_excl])
    #             # print(gas.X[gas.species_index(compound.name)])
    #             # print(gas.X[gas.species_index(compound.name)])
    ############################################################################################
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
    X_CH4_profile.append(reactors[k].thermo.X[gas.species_index('CH4')])
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
ax2.set_ylabel('CH4 profile')

ax3.plot(z_vec, H2_profile)
ax3.set_ylabel('H2 profile')
ax2.set_xlabel('Reactor length')

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))