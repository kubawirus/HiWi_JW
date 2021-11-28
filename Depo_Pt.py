import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import V_N
import Reduction_mech
import json
"""
ZIEL:

1. Einpflegen von Mechanismus 13 in der Gasphase DONE

2. Propandehydrierung mit einfacher Kinetik auf der Oberfläche DONE

3. Testen ob das funktioniert (mit verändern zb. der Konzentrationen und Temperaturen -> sinnvoller Einfluss auf Reaktion/Konzentrationen?) DONE

4. dann: Deaktivierung durch Kohlenstoffablagerungen (Verringerung von der Oberfläche durch große Kohlenwasserstoffspezies > c6h6) ...JETZT

PROBLEMS:

5. Alle Liste als ein Dict zu erschaffen

6. Es scheint so als die Oberfläche keinen Einfluss auf die Reaktionen hat. Nachdem ich die surf pro vol verringert habe
    haben sich die Verläufe von C3H8 usw. nicht geändert. 

7. Oberfläche Einflüss zu untersuchen.
"""

##################### INPUT DATA FROM TERMINAL ########################

# n_steps = int(input('Number of reactors in the pfr system. The number must be divisible by 4!\n'))
# while n_steps % 4 != 0:
#     print("The number must be divisible by 4!\n")
#     n_steps = int(input("Number of reactors:\n"))
# iters = int(input('Number of iterations for the pfr-reactor cylcle?\n'))
# mechanism = int(input('Press 0 for automatic gri_30, 1 to choose reduced gri30 or Press 2 to choose mech_13\n'))
# remove = int(input('Should all heavy components be removed after every reactor? 1- yes 0- no\n'))
n_steps = 40
mechanism = 2
remove = 0
# starts counting time in which program runs
start_time = time.time()

#################### INPUT DATA #######################################
# relative temperature [K]
T_0_n = 273.0
# inlet temperature [K]
T_0 = 473.0
# wall temperature [K]
T_wall = 600 + 273.15
# constant pressure [Pa]
pressure = ct.one_atm
# flow rate [m3/s] volumen flow standarized for 273K
vol_flow_rate_N = 1.6664e-6
# composition
composition_0 = 'C3H8:10, H2:1'
# composition_0 = 'CH4:1, O2:1.5, AR:0.1'
initial_state = T_0, pressure, composition_0
reactive_state = T_wall, pressure, composition_0

# reaction mechanism file.yaml
# reaction_mech = 'mech_13.yaml'

# reaction mechanism for surface reaction
# reaction_mech_surf = 'gri30_PtSurf.yaml'
reaction_mech_surf = 'Propan_surface.yaml'

# Reactor geometry
length = 0.1  # *approximate* PFR length [m]
area = 0.00024  # cross-sectional area of reactor [m**2]
height = 0.006  # [m]
depth = 0.04  # [m]
porosity = 0.6 # [-] It's not a porosity per se, it defines what fraction of reactor will be surface.
cat_area_per_vol = 1000   # Catalyst particle surface area per unit volume [1/m] What is a good value?
area_react = area * (1-porosity)
area_cat = area * porosity
beta = 1 # How much from the main flow (gas reactor) will flow to the surface reactor, 1 = whole mass stream
###################### INITIAL CALCULATIONS & CREATING OBJECTS #############################

# Here the reaction mechanism will be chosen. The default mechanism is gri30.yaml
if mechanism > 0:
    # Calling a function from "Reduction_mech.py" to choose the desired mechanism after reduction.
    # Change of mechanism is avaliable through this named file.

    gas, reaction_mech, gas_surf = Reduction_mech.gas_funct(mechanism, reactive_state)
else:
    reaction_mech = 'gri30.yaml'
    # Import the gas model and set the initial conditions
    gas = ct.Solution(reaction_mech)
    # Import gas model for surface reaction
    gas_surf = ct.Solution(reaction_mech_surf, 'gas')

# # Import gas model
# gas = ct.Solution(reaction_mech)
# State with standard temp 273 for vol. calculation
gas.TPX = T_0_n, pressure, composition_0

gas_surf.TPX = initial_state
# spec_surf = gas_surf.species_names
# spec_gas = gas.species_names

#import the surface model
surf = ct.Interface(reaction_mech_surf, 'Pt_surf', adjacent=[gas_surf])
surf.TP = T_wall, pressure # Which temperature cat should have?
cov = surf.coverages

# calculate a flow rate & velocity
# vol_flow_rate = V_N.vol_norm(vol_flow_rate_N, gas, T_0, pressure)  # [m3/s]
vol_flow_rate = vol_flow_rate_N * T_0_n/T_0
u_0 = vol_flow_rate / area_react  # inflow velocity [m/s]
u_0_surf = vol_flow_rate / area_cat
# Back to initial state
gas.TPX = initial_state
# gas()

# Calculate other data
# mass flow rate
mass_flow_rate_0 = vol_flow_rate * gas.density
# reactor slice length
dz = length / n_steps
# Volume of every slice reactor
r_vol = area_react * dz
# Side wall area of slice reactor (for heat exchange) Q. AT THIS POINT
A_wall = 2 * height * dz + 2 * depth * dz

# cat volume and catalyst surface area in one reactor
cat_vol = area_cat * dz  # This reactor volume is smaller as it represents only near wall zone
cat_surf = cat_vol * cat_area_per_vol

# create a new reactor for gas phase
reactor = ct.IdealGasReactor(contents=gas, energy='on')
reactor.volume = r_vol
# create a reactor for surface reaction
reactor_surf = ct.IdealGasReactor(contents=gas_surf, energy = 'on') #CHANGE:contents=gas_surf
reactor_surf.volume = cat_vol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream = ct.Reservoir(gas, name='upstream')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream = ct.Reservoir(gas, name='downstream')

# Add the reacting surface to the reactor. The area is set to the desired
# catalyst area in the reactor.
rsurf = ct.ReactorSurface(surf, reactor_surf, A=cat_surf)

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
m = ct.MassFlowController(upstream, reactor)
m.mass_flow_rate = mass_flow_rate_0
# Between reactor and surf_reactor
m_surf = ct.MassFlowController(reactor, reactor_surf)
m_surf.mass_flow_rate = mass_flow_rate_0 * beta  # * 0.01 How much from the main flow will flow to the surface

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v = ct.PressureController(reactor, downstream, master=m, K=1e-6)
# Between reactor and surf_reactor
v_surf = ct.PressureController(reactor_surf, reactor, master=m_surf, K=1e-6)

# create reservoirs for heat exchange with gas2 and for enviroment temp gas 3
gas_600 = ct.Solution('liquidvapor.yaml', 'water')
gas_600.TPX = T_wall, 5e+5, 'H2O:1'
# gas2()
heat_reserv = ct.Reservoir(gas_600)

gas_200 = gas_600
gas_200.TPX = 473.0, 5e+5, 'H2O:1'
# gas3()
env_reserv = ct.Reservoir(gas_200)

################## Initializing some lists for holding results  ##################

# One dicts for all results:
result_dict_gas = {}
list_for_results_gas = ["length", "moles", "mflow", "mol_flow", "u", "state_list", "mass"]
for i in list_for_results_gas:
    result_dict_gas[i] = []

result_dict_surf = {}
list_for_results_surf = ["mflow_surf", "u_surf", "state_list_surf", "Pt(S)", "H(S)"]
for i in list_for_results_surf:
    result_dict_surf[i] = []

# Length of the whole reactor
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
# This list is for saving the next reactor state. After every loop it will be reconstructed from this data.
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
# t_react = 0
t_react_sum = 0
t_react_gas = 0
t_react_surf = 0

# Surface Reactor
# mass flows
mflow_surf = [0] * (n_steps + 1)
mflow_surf[0] = mass_flow_rate_0 * beta
# velocity through the reactor
u_surf = [0] * (n_steps + 1)
u_surf[0] = u_0_surf # mflow_surf[0]/(reactor_surf.density * (cat_vol/dz))
# This list is for saving the next reactor state. After every loop it will be reconstructed from this data.
state_list_surf = [0] * (n_steps + 1)
state_list_surf[0] = reactor_surf.thermo.TPY
# pressure in every reactor
p_profile_surf = [0] * (n_steps + 1)
p_profile_surf[0] = reactor_surf.thermo.P
# temp profile
T_profile_surf = [0] * (n_steps + 1)
T_profile_surf[0] = reactor_surf.T


# Create a dictionary for storing deposition components and their moles
depo = {}
depo_r = {}
mass_depo = [0] * (n_steps + 1)
sum_depo = 0.0

# Lists for coverages
cov_dict = {}
# cov_list = [0]*(n_steps+1)
for i, j in enumerate (surf.species_names) :
    cov_dict[j] = []
#     cov_dict[j].append(cov[i])
Pt_cov = [0] * (n_steps+1)
Pt_cov[0] = cov[0]
H_cov = [0] * (n_steps+1)
H_cov[0] = cov[1]

surf_area_n = []
surf_area_n.append(rsurf.area)

prod_rates_surf = {}
prod_rates_gas = {}
##############################################################################

################### SIMULATION & ITERATION #############################

for n in range(n_steps):
    # delete all existing walls and previos state in the reactor
    reactor.walls.clear()
    # create wall for heat exchange. What Wall area [m^2] and Overall heat transfer coefficient [W/m^2]??
    if n < (3 * (n_steps // 4)):
        ct.Wall(heat_reserv, reactor_surf,  A=A_wall, U=10)
    else:
        ct.Wall(env_reserv, reactor_surf,  A=A_wall, U=10)

    # Reduce & update the cat area
    # rsurf = ct.ReactorSurface(surf, reactor_surf, A = surf_area_n[-1])
    # cat_surf = surf_area_n[-1]
    # create a simulation object
    sim = ct.ReactorNet([reactor, reactor_surf])

    # Set the state of the reservoir to match that of the previous reactor
    gas.TPY = state_list[n]
    # gas_surf.TPY = T_wall, state_list[n][1], state_list[n][2]
    gas_surf.TPY = state_list[n]
    # gas_surf.TPY = # state list or dict
    # surf.coverages = # some list for covarages

    upstream.syncState()
    # integrate the reactor forward in time until steady state is reached
    sim.reinitialize()
    sim.advance_to_steady_state()

    t_react_gas += reactor.mass / mflow[n]
    t_react_surf += reactor_surf.mass / mflow_surf[n]

    # temporal variables
    T = reactor.T
    P = reactor.thermo.P
    Y = reactor.Y

    # update list of pressure and temp
    p_profile[n + 1] = P
    T_profile[n + 1] = T
    # result_dict_gas["p_profile"].append(P)
    # result_dict_gas["T_profile"].append(T)

    # update list of pressure and temp
    p_profile_surf[n + 1] = reactor_surf.thermo.P
    T_profile_surf[n + 1] = reactor_surf.thermo.T
    # result_dict_surf["p_profile_surf"].append(reactor_surf.thermo.P)
    # result_dict_surf["T_profile_surf"].append(reactor_surf.thermo.T)

    # update coverages
    for i, j in enumerate(surf.species_names):
         cov_dict[j].append(surf.coverages[i])

    ############## REMOVE ONE COMPONENT FROM THE REACTOR NET ###########

    # remove water after every reactor
    if remove == 1:

        # These are entry data if no compounds is to be removed
        mass[n + 1] = reactor.mass
        mass_flow_now = mflow[n]
        mflow[n + 1] = mass_flow_now
        mol_flow[n + 1] = mass_flow_now / reactor.thermo.mean_molecular_weight
        state_list[n + 1] = (T, P, Y)
        mass_depo[n] = 0.0

        # For every compound in the reactor show only these with mol. weight ...
        # These will be set to 0.0
        for j, compound in enumerate(reactor.thermo.species()):
            # Important condition is stated here. Only components heavier than 80 g/mol are going to be removed
            if reactor.thermo.molecular_weights[j] > 80.00:

                name_Y_excl = compound.name
                Y_excl = gas.species_index(compound.name)
                # If this component doesn't exist on the list, add the name and value 0.
                if name_Y_excl not in depo.keys():
                    depo[name_Y_excl] = 0.0
                    depo_r[name_Y_excl] = []

                # Mass in reactor
                mass[n + 1] -= reactor.mass * Y[Y_excl]
                # Part of the mass flow that is deposited in one slice reactor
                mass_flow_rate_down = mass_flow_now * Y[Y_excl]
                # Sum of deposited mass across reactor
                mass_depo[n] += mass_flow_rate_down
                # And whole summed deposited mass
                sum_depo += mass_flow_rate_down
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
        m.mass_flow_rate = mflow[n+1]
        mflow_surf[n+1] = mflow[n+1] * porosity
        m_surf.mass_flow_rate = mflow_surf[n+1]
        u[n+1] = mflow[n+1]/(reactor.density * area_react)
        u_surf[n+1] = mflow_surf[n+1]/(reactor_surf.density * (cat_vol/dz))
        # state_list_surf[n+1] = (reactor_surf.T, reactor_surf.thermo.P, reactor_surf.Y)
        upstream.syncState()

        ######### INACTIVATION OF REACTIVE SURFACE DUE TO DEPOSITION #############

        ## Problem: How to bind the deacrising activity of surface with increasing deposition of heavy compounds.

        # rsurf.area -= mass_depo[n] * 0.5e5
        # surf_area_n.append(rsurf.area)

        ###########################################################################
    else:
        mass[n + 1] = reactor.mass
        mflow[n + 1] = mflow[n]
        mflow_surf[n+1] = mflow[n+1] # *porosity
        u[n + 1] = mflow[n+1]/(reactor.density * area_react)
        u_surf[n + 1] = mflow_surf[n + 1] /( reactor_surf.density * area_cat)
        mol_flow[n + 1] = mflow[n + 1] / reactor.thermo.mean_molecular_weight
        state_list[n + 1] = (T, P, Y)

        #Baustelle!
        # mflow_temp = result_dict_gas["mflow"][-1]
        # result_dict_gas["mass"].append(reactor.mass)
        # result_dict_gas["mflow"].append(mflow_temp)
        # result_dict_surf["mflow"].append(mflow_temp*0.01)
        # result_dict_gas["u"].append(result_dict_gas["mflow"])

        # state_list_surf[n+1] = (reactor_surf.T, reactor_surf.thermo.P, reactor_surf.Y)
        upstream.syncState()

    # Fill two dictionaries with net production rates for surface and for gas phase.
    len_cov = rsurf.coverages.size
    for i, j in enumerate(rsurf.kinetics.kinetics_species_names):
        if j not in prod_rates_surf.keys():
            prod_rates_surf[j] = []
        prod_rates_surf[j].append(rsurf.kinetics.net_production_rates.T[i])

    for i, j in enumerate(reactor.kinetics.kinetics_species_names):
        if j not in prod_rates_gas.keys():
            prod_rates_gas[j] = []
        prod_rates_gas[j].append(reactor.kinetics.net_production_rates.T[i])
##############################################################################

######################## CREATE SOME RESULTS ###############################

# show gas properties
# gas()

# Create a file to store deposited compounds and their mass flow rate
if remove == 1:
    f = open("Deposition.txt", "w")
    f.write("Comp.: MFR[kg/s]:\n\n")
    for key in depo:
        f.write(key + ": " + str(depo[key]) + "\n")
    f.close()

    # Look for the compound that has a biggest sum of deposition

    # Compound with biggest deposition
    depo_comp_name = max(depo,key=depo.get)
    print("The largest fraction of deposition is ", depo_comp_name)
    # Create a list where the deposited mass of this compound is stored along the reactor
    depo_comp_reactor = []
    depo_comp_reactor.append(0.0)
    for i in range(n_steps):
        depo_comp_reactor.append(depo_r[depo_comp_name][i])


################# Compare Results in matplotlib ###########################
z = z_vec[1:]
y = state_list[1:]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Plots for evaluation of reactor T & p')
ax1.plot(z_vec, p_profile)  # [s[1] for s in state_list]
ax1.set_ylabel('Pressure [Pa]')

ax2.plot(z_vec, T_profile)
ax2.set_ylabel('Temperature [K]')

ax3.plot(z_vec, u)
ax3.set_ylabel('Velocity [m/s]')

ax3.set_xlabel('Reactor length [m]')

fig2, (ax4, ax5, ax6) = plt.subplots(3, 1)
fig2.suptitle('Flow rates & mass across the reactor')
ax4.plot(z_vec, mflow)
ax4.set_ylabel('mass flow rate [kg/s]')

ax5.plot(z_vec, mol_flow)
ax5.set_ylabel('mol flow rate [mol/s]')

ax6.plot(z_vec, mass)
ax6.set_ylabel('mass [kg]')

ax6.set_xlabel('Reactor length [m]')
if remove == 1:
    fig3, (ax7, ax8) = plt.subplots(2,1)
    fig3.suptitle('Deposition across the reactor')

    ax7.plot(z_vec, mass_depo)
    ax7.set_ylabel('sum of depo [kg/s]')

    ax8.plot(z_vec, depo_comp_reactor)
    ax8.set_ylabel(depo_comp_name + ' depo [kg/s]')

    # ax16.plot(z_vec, surf_area_n)
    # ax16.set_ylabel('surface area')

    ax8.set_xlabel('Reactor length [m]')
# Plot for species fractions thru reactor
fig4, (ax9,ax10) = plt.subplots(2,1)
fig4.suptitle('Fractions of species')
# Define list for desired species
c3h8 = []
c3h6 = []

# Fulfill the list
# for i, val in enumerate(state_list):
#     c3h8.append(val[2][gas.species_index('CH4')])
#     c3h6.append(val[2][gas.species_index('H2')])

for i, val in enumerate(state_list):
    c3h8.append(val[2][gas.species_index('C3H8')])
    c3h6.append(val[2][gas.species_index('C3H6')])

# ax9.plot(z_vec, c3h8)
# ax9.set_ylabel('Mass fraction of CH4')
# ax10.plot(z_vec, c3h6)
# ax10.set_ylabel('Mass fraction of H2  ')
# ax10.set_xlabel('Reactor length [m]')

ax9.plot(z_vec, c3h8)
ax9.set_ylabel('Mass fraction of C3H8')
ax10.plot(z_vec, c3h6)
ax10.set_ylabel('Mass fraction of C3H6  ')
ax10.set_xlabel('Reactor length [m]')

# Plots for surface reactor
fig5, (ax11, ax12, ax13) = plt.subplots(3, 1)
fig5.suptitle('Plots for evaluation of surface reactor u & mflow')

ax11.plot(z_vec, T_profile_surf)
ax11.set_ylabel('Temperature [K]')

ax12.plot(z_vec, u_surf)
ax12.set_ylabel('Velocity [m/s]')

ax13.plot(z_vec, mflow_surf)
ax13.set_ylabel('mass flow rate [kg/s]')
ax13.set_xlabel('Reactor length [m]')

# Plot for species fractions thru surface reactor
fig6, (ax14,ax15) = plt.subplots(2,1)
fig6.suptitle('Coverages in surface reactor')
cov1 = cov_dict.get('PT(S)')
cov2 = cov_dict.get('H(S)')
ax14.plot(z, cov1)
ax14.set_ylabel('Cov PT(S)')
ax15.plot(z, cov2)
ax15.set_ylabel('Cov H(S)')
ax15.set_xlabel('Reactor length [m]')

fig7, (ax16, ax17) = plt.subplots(2,1)
fig7.suptitle('Net production rates on the surface')
ax16.plot(z, prod_rates_surf['C3H8'])
ax17.plot(z, prod_rates_surf['C3H6'])
ax16.set_ylabel('C3H8')
ax17.set_ylabel('C3H6')
ax17.set_xlabel('Reactor length [m]')

fig8, (ax18, ax19) = plt.subplots(2,1)
fig8.suptitle('Net production rates in gas phase')
ax18.plot(z, prod_rates_gas['C3H8'])
ax19.plot(z, prod_rates_gas['C3H6'])
ax18.set_ylabel('C3H8')
ax19.set_ylabel('C3H6')
ax19.set_xlabel('Reactor length [m]')

# plt.savefig('Depo_Pt_Plot.png')
plt.show()
# fig4.show()
# fig6.show()
# fig7.show()
# fig8.show()

##################################################################

# Check if sum depo is equal to mass flow rate difference:
print("\n sum depo: ", sum_depo, " mfr difference: ", mass_flow_rate_0 - mflow[n_steps])

t_react_sum = t_react_surf + t_react_gas  # Frage?!
print("\nResidence time in reactive Part = ", t_react_sum, " s")
print("\nResidence time in reactive Part of gas Reactor = ", t_react_gas, " s")
print("\nResidence time in reactive Part of surface reactor = ", t_react_surf, " s")


print("\n--- %s seconds ---" % (time.time() - start_time))
