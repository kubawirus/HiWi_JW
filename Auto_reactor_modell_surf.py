import cantera as ct
# import numpy as np
# import matplotlib.pyplot as plt
"""
In this automated version of reactor model, the number of rectors is free to choose but it must be
divisible by 3. In the first third of reactor there is only heat exchange through the wall, in the second third
the catalytic surface reaction takes place and in the last third comes also to the heat exchange through the 
wall with the first third part of reactor. 

-Without plots
-Outcome (Temp at the end) for 12 reactors is the same as for the "Finale_Version" from Moritz.
"""
#######################################################################
# Input Parameters
#######################################################################

# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.
# IT MUST DIVISBLE BY 3!!!!
n_steps = int(input('Number of reactors in the pfr system. The number must be divisible by 3!\n'))
steps = int(input('Number of iterations for the pfr-reactor cylcle?\n'))
walls = int(input('Press 1 to activate walls for heat transfer\n'))

# unit conversion factors to SI
cm = 0.01
minute = 60.0

T_0 = 1000.0  # inlet temperature [K]
pressure = ct.one_atm  # constant pressure [Pa]
composition_0 = 'CH4:1, O2:1.5, AR:0.1'

#catalyst
# length = 0.3 * cm  # Catalyst bed length
# area = 1.0 * cm**2  # Catalyst bed area
cat_area_per_vol = 1000.0 / cm  # Catalyst particle surface area per unit volume
initial_coverage_catalyst= 'O(S):0.00, PT(S):0.01, H(S):0.99' #Zeile selbst hizugefügt. Wert aus Mechanismus

# Warum solche Angaben??
length = 1.5e-7  # *approximate* PFR length [m]
u_0 = .006  # inflow velocity [m/s]
area = 1.e-4  # cross-sectional area [m**2]

reaction_mechanism1 ='methane_pox_on_pt.cti'

# import the gas model and set the initial conditions
gas = ct.Solution(reaction_mechanism1, 'gas')
gas.TPX = T_0, pressure, composition_0

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
reactors = [ct.IdealGasReactor(gas) for i in range (n_steps)]

# Add a volume to reactors
for i, reactors_obj in enumerate(reactors):
    reactors_obj.volume=r_vol

# Add walls between reactors on the opposite side of system
if walls == 1:
    n_walls = n_steps//3
    # print(n_walls)

    for i in range(n_walls):
        ct.Wall(reactors[i], reactors[n_steps-1-i], A=0.001, U=1.2)
    print('Walls are active!')
else:
    print('Wall are not active!')

# create a reservoir to represent the reactor immediately upstream
upstreams = [ct.Reservoir(gas) for i in range (n_steps)]

# create a downstream reservoir
downstreams = [ct.Reservoir(gas) for i in range (n_steps)]

# Add the reacting surface to the reactor. The area is set to the desired
# catalyst area in the reactor.
rsurfs = []
for i in range (n_steps//3,2*(n_steps//3)):
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

sim = ct.ReactorNet(reactors)

# "Länge" des Gesamtreaktors nach jedem Teilreaktor" ??
z_vec=[]
for i in range(1,n_steps+1):
    z_vec.append(i*dz)

# Zustände in den Reaktoren
dict_states = {}
for i in range(n_steps):
    dict_states[i] = ct.SolutionArray(reactors[i].thermo)

# Create an empty dictionary to store Temperatures of individual reactors
dict_T = {}
# Adding arrays to this dictionary
for i in range(n_steps):
    dict_T[i] = []

#Create Lists for evaluation
temp_profile=[]
p_profile=[]
X_CH4_profile=[]
X_O2_profile=[]
X_CO2_profile=[]
X_H2O_profile=[]

# Create a list for velocity for every reactor
u = []
# For residence time in each reactor
t_r = []
# Whole time
t = []
# Density in each reactor
r_dens= []

# The outer loop iterate through the given number of iterations for every pfr cycle
for n in range(steps):
    # This loop iterate through all reactors
    for i in range(n_steps):

        # Set the state of the reservoir to match that of the previous reactor
        if i > 0:
            gas.TDY = reactors[i-1].thermo.TDY
        else:
            gas.TPX = T_0, pressure, composition_0

        upstreams[i].syncState()

        # integrate the reactor forward in time until steady state is reached
        sim.reinitialize()
        sim.advance_to_steady_state()

        # compute velocity and transform into time
        u.append(mass_flow_rate/area/reactors[i].density)
        t_r.append(reactors[i].mass / mass_flow_rate)  # residence time in this reactor
        t.append(t_r[i])
        r_dens.append(reactors[i].thermo.density)

        # Save temperature & thermo.state in each Reaktor
        for j in range(n_steps):
            dict_T[j].append(reactors[j].T)
            dict_states[j].append(reactors[j].thermo.state)

    # Creating lists with temp, pressure and species conc. for diagrams.
    for k in range(n_steps):
        temp_profile.append(reactors[k].T)
        p_profile.append(reactors[k].thermo.P)
        X_CH4_profile.append(reactors[k].thermo.X[gas.species_index('CH4')])
        X_O2_profile.append(reactors[k].thermo.X[gas.species_index('O2')])
        X_CO2_profile.append(reactors[k].thermo.X[gas.species_index('CO2')])
        X_H2O_profile.append(reactors[k].thermo.X[gas.species_index('H2O')])

    print('Reaktordurchlauf:\t', n+1)

# At the outlet from Reactor
gas.TDY = reactors[n_steps - 1].thermo.TDY
gas()











