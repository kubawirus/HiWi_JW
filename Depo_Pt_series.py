import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import V_N
import Reduction_mech
import json

#-------------------- INPUT DATA FROM TERMINAL --------------------#

# n_steps = int(input('Number of reactors in the pfr system. The number must be divisible by 4!\n'))
# while n_steps % 4 != 0:
#     print("The number must be divisible by 4!\n")
#     n_steps = int(input("Number of reactors:\n"))
# iters = int(input('Number of iterations for the pfr-reactor cylcle?\n'))
# mechanism = int(input('Press 0 for automatic gri_30, 1 to choose reduced gri30 or Press 2 to choose mech_13\n'))
# inactiv = int(input('Should surface be decreased due to deposition? 1- yes'))
n_steps = 40
mechanism = 2
remove = 1
inactiv = 1 # inactiv = 1 means you want deposition to have an influence on surface of react_surf
# starts counting time in which program runs
start_time = time.time()

#-------------------- INPUT DATA --------------------#
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
# definition of initial and reactive state
initial_state = T_0, pressure, composition_0
reactive_state = T_wall, pressure, composition_0
# Surface deactivation coefficient - is used to artificially enlarge the influence of deposition on a surface
alpha = 10e5

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
porosity = 0.6 # [-] It's not a porosity per se, it defines what fraction of reactor will be surface reactor.
cat_area_per_vol = 1000   # Catalyst particle surface area per unit volume [1/m] What is a good value?
area_react = area * (1-porosity)
area_cat = area * porosity
beta = 1 # How much from the main flow (gas reactor) will flow to the surface reactor, 1 = whole mass stream

#-------------------- INITIAL CALCULATIONS & CREATING OBJECTS --------------------#

# Here the reaction mechanism will be chosen. The default mechanism is gri30.yaml.
# The mechanism will be reduced, so that only most active reaction and reactants are present in the system.
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

# SSet gas into initial state
gas.TPX = initial_state
# State with initial temp
gas_surf.TPX = initial_state

# import the surface model
surf = ct.Interface(reaction_mech_surf, 'Pt_surf', adjacent=[gas_surf])
surf.TP = T_0, pressure
cov = surf.coverages

# calculate a flow rate & velocity
# vol_flow_rate = V_N.vol_norm(vol_flow_rate_N, gas, T_0, pressure)  # [m3/s]
vol_flow_rate = vol_flow_rate_N * T_0_n/T_0
u_0 = vol_flow_rate / area_react  # inflow velocity [m/s]
u_0_surf = vol_flow_rate / area_cat

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
m_surf.mass_flow_rate = mass_flow_rate_0 * beta  # beta - How much from the main flow will flow to the surface

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

#-------------------- Initializing dictionaries for holding results --------------------#

# One dicts for all results:
result_dict_gas = {}
list_for_results_gas = ["length", "moles", "mflow", "mol_flow", "u", "state_list", "mass", "pressure", "temp"]
for i in list_for_results_gas:
    result_dict_gas[i] = []

result_dict_surf = {}
list_for_results_surf = ["mflow_surf", "u_surf", "state_list_surf", "surf_area_n"]
for i in list_for_results_surf:
    result_dict_surf[i] = []

#---- Coverages ----#
len_cov = rsurf.coverages.size
result_dict_surf ["cov"] = {}
for i, j in enumerate (surf.species_names) :
    result_dict_surf["cov"][j] = []

# ---- Deposition ----#
# Create a dictionary for storing deposition components and their moles
result_dict_surf ["depo"] = {}
result_dict_surf ["depo_r"] = {}
result_dict_surf ["mass_depo"] = []
result_dict_surf ["sum_depo"] = 0.0

result_dict_surf["mass_depo"].append(0)

# ---- Net Production Rates ----#
result_dict_surf["prod_rates_surf"] = {}
result_dict_gas["prod_rates_gas"] = {}


# Length of the whole reactor
for i in range(n_steps + 1):
    result_dict_gas["length"].append(i * dz)

# Time in reactive part
t_react_sum = 0
t_react_gas = 0
t_react_surf = 0

#Update all lists for gas reactor

#---- Gas Reactor ----#
result_dict_gas["state_list"].append(reactor.thermo.TPY)
result_dict_gas["pressure"].append(reactor.thermo.P)
result_dict_gas["temp"].append(reactor.T)
result_dict_gas["mflow"].append(mass_flow_rate_0)
result_dict_gas["u"].append(u_0)
result_dict_gas["mass"].append(reactor.mass)

#---- Surface Reactor ----#
result_dict_surf["mflow_surf"].append(mass_flow_rate_0 * beta)
result_dict_surf["u_surf"].append(u_0_surf)
result_dict_surf["surf_area_n"].append(cat_surf)

#---- Time ----#
t_react_gas += reactor.mass / mass_flow_rate_0
t_react_surf += reactor_surf.mass / mass_flow_rate_0 * beta

# Initialize state list of reactor
STATE_LIST = reactor.thermo.TPY

#----------------------------------------------------------------#

#-------------------- SIMULATION & ITERATION --------------------#
for n in range(n_steps):
    # delete all existing walls and previos state in the reactor
    reactor.walls.clear()
    # create wall for heat exchange. What Wall area [m^2] and Overall heat transfer coefficient [W/m^2]??
    if n < (3 * (n_steps // 4)):
        ct.Wall(heat_reserv, reactor_surf,  A=A_wall, U=10)
    else:
        ct.Wall(env_reserv, reactor_surf,  A=A_wall, U=10)

    # create a simulation object
    sim = ct.ReactorNet([reactor, reactor_surf])

    # Set the state of the reservoir to match that of the previous reactor
    gas.TPY = STATE_LIST
    gas_surf = STATE_LIST

    upstream.syncState()

    # integrate the reactor forward in time until steady state is reached
    sim.reinitialize()
    sim.advance_to_steady_state()

    # temporal variables, these are mit capital letters
    T = reactor.T
    P = reactor.thermo.P
    Y = reactor.Y
    M = reactor.mass
    M_FLOW = m.mass_flow_rate
    STATE_LIST = (T, P, Y)

    # update list of pressure and temp
    result_dict_gas["pressure"].append(P)
    result_dict_gas["temp"].append(T)

    # update coverages
    for i, j in enumerate(surf.species_names):
         result_dict_surf["cov"][j].append(surf.coverages[i])

    #-------------------- REMOVE ONE COMPONENT FROM THE REACTOR NET --------------------#

    # remove water after every reactor
    if remove == 1:

        # These are entry data if no compounds is to be removed
        result_dict_gas["mass"].append(M)
        result_dict_gas["mflow"].append(M_FLOW)
        # mol_flow[n + 1] = mass_flow_now / reactor.thermo.mean_molecular_weight
        result_dict_gas["state_list"].append((T,P,Y))
        result_dict_surf["mass_depo"].append(0.0)

        # For every compound in the reactor show only these with mol. weight ...
        # These will be set to 0.0
        for j, compound in enumerate(reactor.thermo.species()):
            # Important condition is stated here. Only components heavier than 80 g/mol are going to be removed
            if reactor.thermo.molecular_weights[j] > 80.00:

                name_Y_excl = compound.name
                Y_excl = gas.species_index(name_Y_excl)
                # If this component doesn't exist on the list, add the name and value 0.
                if name_Y_excl not in result_dict_surf["depo"].keys():
                    result_dict_surf["depo"][name_Y_excl] = 0.0
                    result_dict_surf["depo_r"][name_Y_excl] = []

                # Mass in reactor
                M -= reactor.mass * Y[Y_excl]
                # Part of the mass flow that is deposited in one slice reactor
                mass_flow_rate_down = M_FLOW * Y[Y_excl]
                # Sum of deposited mass across reactor
                result_dict_surf["mass_depo"][-1] += mass_flow_rate_down
                # And whole summed deposited mass
                result_dict_surf["sum_depo"] += mass_flow_rate_down
                # Set new mass flow rate in upstream
                M_FLOW -= mass_flow_rate_down
                # Set new mol flow rate in the upstream
                # mol_flow[n + 1] -= mass_flow_rate_down / gas.molecular_weights[Y_excl]

                # Set new mass fraction after deposition (later it can be a function)
                Y[Y_excl] = 0.0
                STATE_LIST = (T, P, Y)


                # Save deposited compound and moles in a dict.
                result_dict_surf["depo"][name_Y_excl] += mass_flow_rate_down
                result_dict_surf["depo_r"][name_Y_excl].append(mass_flow_rate_down)

        # Set a new value for mass flow in the upstream
        m.mass_flow_rate = M_FLOW
        m_surf.mass_flow_rate = M_FLOW * beta
        result_dict_gas["mass"][-1] = M
        result_dict_gas["mflow"][-1] = M_FLOW
        result_dict_surf["mflow_surf"][-1] = M_FLOW * beta
        result_dict_gas["u"].append(M_FLOW/(reactor.density * area_react))
        result_dict_surf["u_surf"].append( (M_FLOW * beta )/ (reactor.density * area_cat))
        result_dict_gas["state_list"][-1] = ((T, P, Y))
        upstream.syncState()

        #-------------- DEACTIVATION OF REACTIVE SURFACE DUE TO DEPOSITION ---------#
        if inactiv == 1:
        ## Problem: How to bind the deacrising activity of surface with increasing deposition of heavy compounds.
            cat_area_per_vol_new = cat_area_per_vol * max(10e-8, (1 - (result_dict_surf["sum_depo"] * alpha)))
            cat_surf = cat_vol * cat_area_per_vol_new
            rsurf.area = cat_surf
            result_dict_surf["surf_area_n"].append(rsurf.area)

        #----------------------------------------------------------------------------#
    else:
        #Baustelle!
        result_dict_gas["mass"].append(M)
        result_dict_gas["mflow"].append(M_FLOW)
        result_dict_surf["mflow_surf"].append(M_FLOW * beta)
        result_dict_gas["u"].append(M_FLOW/(reactor.density * area_react))
        result_dict_surf["u_surf"].append( (M_FLOW * beta )/ (reactor.density * area_cat))
        result_dict_gas["state_list"].append(STATE_LIST)
        upstream.syncState()

    # Update residence time
    t_react_gas += reactor.mass / M_FLOW
    t_react_surf += reactor_surf.mass / M_FLOW * beta

    # Fill two dictionaries with net production rates for surface and for gas phase.
    for i, j in enumerate(rsurf.kinetics.kinetics_species_names):
        if j not in result_dict_surf["prod_rates_surf"].keys():
            result_dict_surf["prod_rates_surf"][j] = []
        result_dict_surf["prod_rates_surf"][j].append(rsurf.kinetics.net_production_rates.T[i])

    for i, j in enumerate(reactor.kinetics.kinetics_species_names):
        if j not in result_dict_gas["prod_rates_gas"].keys():
            result_dict_gas["prod_rates_gas"][j] = []
        result_dict_gas["prod_rates_gas"][j].append(reactor.kinetics.net_production_rates.T[i])
#----------------------------------------------------------------------------#

#-------------- CREATE SOME RESULTS --------------#

# show gas properties
# gas()

# Create a file to store deposited compounds and their mass flow rate
if remove == 1:
    f = open("Deposition.txt", "w")
    f.write("Comp.: MFR[kg/s]:\n\n")
    for key in result_dict_surf["depo"]:
        f.write(key + ": " + str(result_dict_surf["depo"][key]) + "\n")
    f.close()

    # Look for the compound that has a biggest sum of deposition

    # Compound with biggest deposition
    depo_comp_name = max(result_dict_surf["depo"],key=result_dict_surf["depo"].get)
    print("The largest fraction of deposition is ", depo_comp_name)
    # Create a list where the deposited mass of this compound is stored along the reactor
    depo_comp_reactor = []
    depo_comp_reactor.append(0.0)
    for i in range(n_steps):
        depo_comp_reactor.append(result_dict_surf["depo_r"][depo_comp_name][i])


#-------------- Compare Results in matplotlib --------------#
z = result_dict_gas["length"]
# y = state_list[1:]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Plots for evaluation of reactor T & p')
ax1.plot(z, result_dict_gas["pressure"])  # [s[1] for s in state_list]
ax1.set_ylabel('Pressure [Pa]')

ax2.plot(z, result_dict_gas["temp"])
ax2.set_ylabel('Temperature [K]')

ax3.plot(z, result_dict_gas["u"])
ax3.set_ylabel('Velocity [m/s]')

ax3.set_xlabel('Reactor length [m]')

fig2, (ax4, ax6) = plt.subplots(2, 1)
fig2.suptitle('Flow rates & mass across the reactor')
ax4.plot(z, result_dict_gas["mflow"])
ax4.set_ylabel('mass flow rate [kg/s]')

ax6.plot(z, result_dict_gas["mass"])
ax6.set_ylabel('mass [kg]')

ax6.set_xlabel('Reactor length [m]')

if remove == 1:
    fig3, (ax7, ax8) = plt.subplots(2,1)
    fig3.suptitle('Deposition across the reactor')

    ax7.plot(z, result_dict_surf["mass_depo"])
    ax7.set_ylabel('sum of depo [kg/s]')

    ax8.plot(z, depo_comp_reactor)
    ax8.set_ylabel(depo_comp_name + ' depo [kg/s]')

    ax8.set_xlabel('Reactor length [m]')

# Plot for species fractions thru reactor
fig4, (ax9,ax10) = plt.subplots(2,1)
fig4.suptitle('Fractions of species')
# Define list for desired species
list_species1 = []
list_species2 = []

# Fulfill the list
for i, val in enumerate(result_dict_gas["state_list"]):
    list_species1.append(val[2][gas.species_index('C3H8')])
    list_species2.append(val[2][gas.species_index('C3H6')])

ax9.plot(z, list_species1)
ax9.set_ylabel('Mass fraction of C3H8',)
ax10.plot(z, list_species2)
ax10.set_ylabel('Mass fraction of C3H6')
ax10.set_xlabel('Reactor length [m]')

# Plots for surface reactor
# fig5, (ax11, ax12, ax13) = plt.subplots(3, 1)
# fig5.suptitle('Plots for evaluation of surface reactor u & mflow')
#
# ax11.plot(z_vec, T_profile_surf)
# ax11.set_ylabel('Temperature [K]')
#
# ax12.plot(z_vec, u_surf)
# ax12.set_ylabel('Velocity [m/s]')
#
# ax13.plot(z_vec, mflow_surf)
# ax13.set_ylabel('mass flow rate [kg/s]')
# ax13.set_xlabel('Reactor length [m]')

# Plot for species fractions thru surface reactor
fig6, (ax14,ax15) = plt.subplots(2,1)
fig6.suptitle('Coverages in surface reactor')
cov1 = result_dict_surf["cov"].get('PT(S)')
cov2 = result_dict_surf["cov"].get('H(S)')
ax14.plot(z[1:], cov1)
ax14.set_ylabel('Cov PT(S)')
ax15.plot(z[1:], cov2)
ax15.set_ylabel('Cov H(S)')
ax15.set_xlabel('Reactor length [m]')

fig7, (ax16, ax17) = plt.subplots(2,1)
fig7.suptitle('Net production rates on the surface')
ax16.plot(z[1:], result_dict_surf["prod_rates_surf"]['C3H8'])
ax17.plot(z[1:], result_dict_surf["prod_rates_surf"]['C3H6'])
ax16.set_ylabel('C3H8')
ax17.set_ylabel('C3H6')
ax17.set_xlabel('Reactor length [m]')

fig8, (ax18, ax19) = plt.subplots(2,1)
fig8.suptitle('Net production rates in gas phase')
ax18.plot(z[1:], result_dict_gas["prod_rates_gas"]['C3H8'])
ax19.plot(z[1:], result_dict_gas["prod_rates_gas"]['C3H6'])
ax18.set_ylabel('C3H8')
ax19.set_ylabel('C3H6')
ax19.set_xlabel('Reactor length [m]')

# plt.savefig('Depo_Pt_Plot.png')
plt.show()
# fig4.show()
# if remove == 1:
#     fig3.show()
# fig7.show()
# fig8.show()

#----------------------------------------------------------------------#
# Linear extrapolation of deposition BAUSTELLE!
'''
Annahmen: 
1. Es findet keine Deaktivierung der Oeberflaeche statt.
2. Es folgt, dass jedes mall genauso so viel Ablagerung produziert wird. 
3. 
'''
t_react_sum = t_react_surf + t_react_gas
mass_depo_ext = ( result_dict_surf["sum_depo"] / t_react_sum ) * 3600 * 1000
print('\nDeposition in reactor ', mass_depo_ext , 'g/h')
# Check if sum depo is equal to mass flow rate difference:
print("\n check sum : ", (result_dict_surf["sum_depo"] - ( mass_flow_rate_0 - result_dict_gas["mflow"][-1])))

print("\nResidence time in reactive Part = ", t_react_sum, " s")
print("\nResidence time in reactive Part of gas Reactor = ", t_react_gas, " s")
print("\nResidence time in reactive Part of surface reactor = ", t_react_surf, " s")


print("\n--- %s seconds ---" % (time.time() - start_time))
