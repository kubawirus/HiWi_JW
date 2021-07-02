
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


#######################################################################
# Input Parameters
#######################################################################
remove = int(input('Should water be removed after every reactor? 1- yes 0- no\n'))

T_0 = 1500.0  # inlet temperature [K]
pressure = ct.one_atm  # constant pressure [Pa]
composition_0 = 'H2:2, O2:1, AR:0.1'
length = 1.5e-7  # *approximate* PFR length [m]
u_0 = .006  # inflow velocity [m/s]
area = 1.e-4  # cross-sectional area [m**2]

# input file containing the reaction mechanism
reaction_mechanism = 'h2o2.yaml'

# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.
n_steps = 100

# import the gas model and set the initial conditions
gas2 = ct.Solution(reaction_mechanism)
gas2.TPX = T_0, pressure, composition_0
mass_flow_rate2 = u_0 * gas2.density * area
dz = length / n_steps
r_vol = area * dz
print(ct.Quantity(gas2).moles)

# create a new reactor
r2 = ct.IdealGasReactor(gas2)
r2.volume = r_vol

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream = ct.Reservoir(gas2, name='upstream')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream = ct.Reservoir(gas2, name='downstream')

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
m = ct.MassFlowController(upstream, r2, mdot=mass_flow_rate2)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v = ct.PressureController(r2, downstream, master=m, K=1e-5)

sim2 = ct.ReactorNet([r2])

# "Länge" des Gesamtreaktors nach jedem Teilreaktor" ??
z_vec = []
for i in range(n_steps + 1):
    z_vec.append(i * dz)
# moles
moles = []
# mass flows
mflow = []
mflow.append(mass_flow_rate2)
# mole flow
mol_flow = []
mol_flow.append(mass_flow_rate2 / r2.thermo.mean_molecular_weight)
# This dictionary is for saving the next reactor state. After every loop it will be reconstructed from this data.
state_list = []
states2 = r2.thermo.TPY
state_list.append(states2)
# mass in every reactor
mass = []
mass.append(r2.mass)
# pressure in every reactor
p_profile = []
p_profile.append(r2.thermo.P)
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
 
8) Połączyć przepływ masowy z ciśnieniem w reaktorze!!
"""
for n in range(n_steps):
    # Set the state of the reservoir to match that of the previous reactor
    gas2.TDY = state_list[n]
    upstream.syncState()
    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()
    sim2.advance_to_steady_state()

    p_profile.append(r2.thermo.P)
    ############## REMOVE ONE COMPONENT FROM THE REACTOR NET ###########3
    if remove == 1:
        # remove water after every reactor
        h20 = gas2.species_index('H2O')
        T, D, Y, = r2.thermo.TDY

        mass.append(r2.mass)
        mol_flow.append(m.mass_flow_rate / r2.thermo.mean_molecular_weight)

        mass_flow_rate_down = m.mass_flow_rate * r2.thermo.Y[h20]
        m.mass_flow_rate = (mflow[n] - mass_flow_rate_down)
        mflow.append(m.mass_flow_rate)

        Y[h20] = 0.0
        state_list.append((T, D, Y))
    else:
        state_list.append(r2.thermo.TDY)
        mass.append(r2.mass)
    ########################################################################


###################################################################
# Compare Results in matplotlib
###################################################################
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Plots for evaluation of reactor')
ax1.plot(z_vec,p_profile) #[s[1] for s in state_list]
ax1.set_ylabel('Pressure')

ax2.plot(z_vec, mflow)
ax2.set_ylabel('mass flow rate')

ax3.plot(z_vec, mol_flow)
ax3.set_ylabel('mol flow rate')
ax2.set_xlabel('Reactor length')

plt.show()
