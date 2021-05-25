"""
A simplistic approach to mechanism reduction which demonstrates Cantera's
features for dynamically manipulating chemical mechanisms.

Here, we use the full GRI 3.0 mechanism to simulate adiabatic, constant pressure
ignition of a lean methane/air mixture. We track the maximum reaction rates for
each reaction to determine which reactions are the most important, according to
a simple metric based on the relative net reaction rate.

We then create a sequence of reduced mechanisms including only the top reactions
and the associated species, and run the simulations again with these mechanisms
to see whether the reduced mechanisms with a certain number of species are able
to adequately simulate the ignition delay problem.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

"""
Aufgabe:
Aus dem mech_13.yaml müssen die Reaktionen ausgeschlossen werden, die kleine (<10^-20) 
Reaktionsrate haben. 
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

gas = ct.Solution('mech_13.yaml')
# Here is the initial conditions are defined
initial_state = 873.0, ct.one_atm, 'C3H8:10, H2:1'

# Run a simulation with the full mechanism
gas.TPX = initial_state
r = ct.IdealGasConstPressureReactor(gas)
sim = ct.ReactorNet([r])

tt = []
TT = []
t = 0.0
# Rmax is the maximum relative reaction rate at any timestep. An Array of size n_reactions full of zeros is created
Rmax = np.zeros(gas.n_reactions)
while t < 0.5:
    t = sim.step()
    tt.append(1000 * t)
    TT.append(r.T)
    rnet = abs(gas.net_rates_of_progress)
    rnet /= max(rnet)
    Rmax = np.maximum(Rmax, rnet)

plt.plot(tt, TT, label='K=All_species, R=All_reac.', color='k', lw=2)

# Get the reaction objects, and sort them so the most active reactions are first
R = sorted(zip(Rmax, gas.reactions()), key=lambda x: -x[0])

# Most active reactions rate
kmax = 10.0**(-20)
# Create an empty list for reactions
reactions = []

# Iterating through every reaction (these are sorted that's why break after else)
# and writing to new list all these reactions that have k > 10^-20
for i, j in enumerate(R):
    if j[0] > kmax:
        reactions.append(j[1])
    else:
        break

print(reactions)

# find the species involved in these reactions. At a minimum, include all
# species in the reactant mixture`
species_names = {'C3H8' , 'H2'}
for reaction in reactions:
    species_names.update(reaction.reactants)
    species_names.update(reaction.products)

# Get the species objects
species = [gas.species(name) for name in species_names]
print(species)

# create the new reduced mechanism
gas2 = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                   species=species, reactions=reactions)
# show gas2
gas2()


# Re-run the ignition problem with the reduced mechanism
gas2.TPX = initial_state
r = ct.IdealGasConstPressureReactor(gas2)
sim = ct.ReactorNet([r])

t = 0.0

tt = []
TT = []
while t < 0.5:
    t = sim.step()
    tt.append(1000 * t)
    TT.append(r.T)

# Checking on the plot if the mechanism looks same with reduced number of reactions.
plt.plot(tt, TT, 'r--', label='K= Reduced spec., R= Reduced react.', lw = 1 )
plt.xlabel('Time (ms)')
plt.ylabel('Temperature (K)')
plt.legend(loc='upper left')
plt.title('Reduced mechanism ignition delay times\n'
          'K: number of species; R: number of reactions')
plt.tight_layout()

plt.show()