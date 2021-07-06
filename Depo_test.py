
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time
import V_N
import Reduction_mech

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
vol_flow_rate_N = 1.6664e-6 # [m3/s] volumen flow standarized for 273K
vol_flow_rate = V_N.vol_norm(vol_flow_rate_N, gas, T_0, pressure) #[m3/s]
length = 0.1  # *approximate* PFR length [m]
area = 0.00024  # cross-sectional area [m**2]
u_0 = vol_flow_rate/area # inflow velocity [m/s]
height = 0.006 #[m]

# Back to initial state
gas.TPX = initial_state
gas()

# Calculate other data
mass_flow_rate_0 = u_0 * gas.density * area
dz = length / n_steps
r_vol = area * dz

# create a new reactor
reactor = ct.IdealGasReactor(contents = gas, energy = 'on')
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
moles = [0]*(n_steps+1)
# mass flows
mflow = [0]*(n_steps+1)
mflow [0] = mass_flow_rate_0
# mole flow
mol_flow = [0]*(n_steps+1)
mol_flow[0] = mass_flow_rate_0 / reactor.thermo.mean_molecular_weight
# This dictionary is for saving the next reactor state. After every loop it will be reconstructed from this data.
state_list = [0]*(n_steps+1)
state_list[0] = reactor.thermo.TDY
# mass in every reactor
mass = [0]*(n_steps+1)
mass[0] = reactor.mass
# pressure in every reactor
p_profile = [0]*(n_steps+1)
p_profile[0] = reactor.thermo.P
# temp profile
T_profile = [0]*(n_steps+1)
T_profile[0] = reactor.T
# density profile
D_profile = [0]*(n_steps+1)
D_profile[0] = reactor.density

# Create a dictionary for storing deposition components and their moles
depo = {}
mass_depo = []

# iterate through the PFR cells

"""
Pomysł jak roziwązać problem zmiennego natężenia przepływu:
1) Stworzyć struktury danych do przechowywania każdego stanu poszczególnego reaktora. 
2) Po każdym reaktorze (n_step) zapisać na końcu nowy stan reaktora oraz nowy strumień masy.
3) Termodynamika gazu oraz przepływ masowy każdego kolejnego reaktora 
    będą definiowane przez poprzedni stan zapisany w odpowiedniej strukturze.
4) Symulacja.
5) Usunięcie związku. Obliczenie nowego przepływu masowego w oparciu o poprzedni.
6) Nowa iteracja może nadpisywać poprzednią. (Może??) 
7) Stworzyć wykresy przepływu masowego przez sieć reaktorów.

Problemy: 
8) Połączyć przepływ masowy z ciśnieniem w reaktorze!! Wydaje się, że ciśnienie jest utrzymywane stałe w reaktorze. 
9) Co ma być stałe, przepływ masowy czy np przepływ obj czyli tylko prędkość?!(Kluczowe!!!)
10) Porównać dwie metody symulacji bez depozycji. 
11) Napisać tak depozycje, żeby kolejne elementy usuwały się "addytywnie"
12) Jak przy depozycji obliczyć gęstość?
13) Stworzyć dict który zapisuje i najlepiej dodaje depozycje z kazdego reaktora!!!
14) Sprawdzic jakie elementy najciezsze mozna deponowac w gri30 zeby bylo szybciej niz w mech_13
"""
# for i in range (iters) :
#     # To reset the old state
#     reactor = ct.IdealGasReactor(contents = gas, energy = 'on')
for n in range(n_steps):
    # delete all existing walls and previos state in the reactor
    # reactor.thermo.TDY.clear()
    reactor.walls.clear()
    # create wall for heat exchange. What Wall area [m^2] and Overall heat transfer coefficient [W/m^2]??
    if (n_steps // 4) <= n < (3 * (n_steps // 4)):
        ct.Wall(heat_reserv, reactor, A=dz * height, U=100)
    else:
        ct.Wall(env_reserv, reactor, A=dz * height, U=100)

    sim2 = ct.ReactorNet([reactor])

    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = state_list[n]
    upstream.syncState()
    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()
    sim2.advance_to_steady_state()

    T = reactor.T
    D = reactor.density
    Y = reactor.Y

    p_profile[n+1] = reactor.thermo.P
    T_profile[n+1] = T
    D_profile[n+1] = D
    ############## REMOVE ONE COMPONENT FROM THE REACTOR NET ###########

    # remove water after every reactor
    if remove == 1:
        # Open txt file to save deposited components and moles
        # f = open("Deposition_DepoTest.txt", "w")
        # f.write("This is a text file where all deposited compounds are saved\n\n")

        mass[n+1] = reactor.mass
        mass_flow_now = D * u_0 * area
        mflow [n+1] = mass_flow_now
        mol_flow[n+1] = mass_flow_now / reactor.thermo.mean_molecular_weight
        state_list[n+1] = (T, D, Y)

        # For every compound in the reactor show only these with mol. weight ...
        # These will be set to 0.0
        for j, compound in enumerate(reactor.thermo.species()):
            if reactor.thermo.molecular_weights[j] > 80.0:

                name_Y_excl = compound.name
                Y_excl = gas.species_index(compound.name)

                if n == 0:
                    dict[str(name_Y_excl)] = 0.0

                mass[n+1] -= reactor.mass * Y[Y_excl]
                # Part of the mass flow that is deposited in one slicereactor
                mass_flow_rate_down = mass_flow_now * Y[Y_excl]
                # # Sum of deposited mass
                # mass_depo[Y_excl] += mass_flow_rate_down
                # Set new mass flow rate in upstream
                mflow[n+1] -= mass_flow_rate_down
                # Set new mol flow rate in the upstream
                mol_flow[n+1] -= mass_flow_rate_down / gas.molecular_weights[Y_excl]

                # Set new mass fraction after deposition (later it can be a function)
                Y[Y_excl] = 0.0
                state_list[n+1] = (T, D, Y)

                # Save deposited compound and moles in a dict.
                depo[str(name_Y_excl)] += mass_flow_rate_down

        # # write these data also to a txt file
        # f.write(name_Y_excl + " " + str(mass_flow_rate_down) + "\n")
        # f.close()
        # Set a new value for mass flow in the upstream
        m.mass_flow_rate = mflow[n+1]
        upstream.syncState()
    else:
        mass [n+1] = reactor.mass
        m.mass_flow_rate = D * u_0 * area
        mflow [n+1] = m.mass_flow_rate
        mol_flow[n+1] = m.mass_flow_rate / reactor.thermo.mean_molecular_weight
        state_list[n+1] = (T, D, Y)
        upstream.syncState()
########################################################################
gas()

###################################################################
# Compare Results in matplotlib
###################################################################
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Plots for evaluation of reactor T & p')
ax1.plot(z_vec,p_profile) #[s[1] for s in state_list]
ax1.set_ylabel('Pressure')

ax2.plot(z_vec, T_profile)
ax2.set_ylabel('Temperature')

ax3.plot(z_vec, D_profile)
ax3.set_ylabel('Density')

ax2.set_xlabel('Reactor length')


fig2, (ax3, ax4, ax5) = plt.subplots(3, 1)
fig2.suptitle('Flow rates & mass across the reactor')
ax3.plot(z_vec, mflow)
ax3.set_ylabel('mass flow rate')

ax4.plot(z_vec, mol_flow)
ax4.set_ylabel('mol flow rate')

ax5.plot(z_vec, mass)
ax5.set_ylabel('mass across reactor')

ax5.set_xlabel('Reactor length')

plt.show()

print("--- %s seconds ---" % (time.time() - start_time))