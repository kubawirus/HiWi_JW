
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
n_steps = 10
iters = 10
#####################################################################

# import the gas model and set the initial conditions
gas = ct.Solution(reaction_mechanism)
gas.TPX = T_0, pressure, composition_0
mass_flow_rate_0 = u_0 * gas.density * area
dz = length / n_steps
r_vol = area * dz
gas()

# Create ideal gas reactors.
reactors = [ct.IdealGasReactor(contents = gas, energy = 'on') for i in range(n_steps + 1)]

# Add a volume to reactors
for i, reactors_obj in enumerate(reactors):
    reactors_obj.volume = r_vol

# create a reservoir to represent the reactor immediately upstream
upstreams = [ct.Reservoir(gas) for i in range(n_steps + 1)]

# create a downstream reservoir
downstreams = [ct.Reservoir(gas) for i in range(n_steps)]


# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
mflows = []
for i in range(n_steps + 1):
    mflows.append(ct.MassFlowController(upstreams[i], reactors[i]))
    mflows[i].mass_flow_rate = mass_flow_rate_0

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference. What K?!
vpress = []
for i in range(n_steps):
    vpress.append(ct.PressureController(reactors[i], downstreams[i], master=mflows[i], K=1e-5))

# "Länge" des Gesamtreaktors nach jedem Teilreaktor" ??
z_vec = []
for i in range(n_steps):
    z_vec.append(i * dz)

# Create a list for velocity for every reactor
p_profile = []
# mass flow rate
mf_rate = []
# moles in reactor
mol_in_react =[]

# iterate through the PFR cells
# The outer loop iterate through the given number of iterations for every pfr cycle
for n in range(iters):
    # This loop iterate through all reactors
    for i in range(n_steps):

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

        if remove == 1:
            # remove water after every reactor
            h20 = gas.species_index('H2O')
            T, D, Y, = downstreams[i].thermo.TDY
            Y[h20] = 0.0

            upstreams[i + 1].thermo.TDY = T, D, Y
            upstreams[i + 1].syncState()

            mol_flow_rate_down = (mflows[i].mass_flow_rate/reactors[i].thermo.mean_molecular_weight) * reactors[i].thermo.Y[h20]
            mflows[i+1].mass_flow_rate -= mol_flow_rate_down


        else:
            upstreams[i + 1].thermo.TPX = downstreams[i].thermo.TPX
            upstreams[i + 1].syncState()



gas.TPX = reactors[n_steps-1].thermo.TPX
gas()
print("moles in 1st reactor: ", reactors[0].mass/reactors[0].thermo.mean_molecular_weight)
print("\nmoles in last reactor: ",reactors[n_steps-1].mass/reactors[n_steps-1].thermo.mean_molecular_weight)

for k in range(n_steps):
    # compute velocity and transform into time
    p_profile.append(reactors[k].thermo.P)
    mf_rate.append(mflows[k].mass_flow_rate)
    mol_in_react.append(mflows[k].mass_flow_rate / reactors[k].thermo.mean_molecular_weight)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
fig.suptitle('Plots for evaluation of reactor')
ax1.plot(z_vec, p_profile)
ax1.set_ylabel('Pressure')

ax2.plot(z_vec, mf_rate)
ax2.set_ylabel('mass flow rate')

ax3.plot(z_vec, mol_in_react)
ax3.set_ylabel('mol flow rate')
ax2.set_xlabel('Reactor length')

plt.show()

#####################################################################