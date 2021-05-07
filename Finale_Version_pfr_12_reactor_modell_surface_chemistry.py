# -*- coding: utf-8 -*-

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

#######################################################################
# Input Parameters
#######################################################################
walls=int(input('Press 1 to activate walls for heat transfer\n'))
steps=int(input('Number of iterations for the pfr-reactor cylcle?\n'))

# unit conversion factors to SI
cm = 0.01
minute = 60.0

T_0 = 1000.0  # inlet temperature [K]
pressure = ct.one_atm  # constant pressure [Pa]
composition_0 = 'CH4:1, O2:1.5, AR:0.1' 

#catalyst
length = 0.3 * cm  # Catalyst bed length
area = 1.0 * cm**2  # Catalyst bed area 
cat_area_per_vol = 1000.0 / cm  # Catalyst particle surface area per unit volume
initial_coverage_catalyst= 'O(S):0.00, PT(S):0.01, H(S):0.99' #Zeile selbst hizugefügt. Wert aus Mechanismus        

length = 1.5e-7  # *approximate* PFR length [m]
u_0 = .006  # inflow velocity [m/s]
area = 1.e-4  # cross-sectional area [m**2]

# input file containing the reaction mechanism
# grundsätzlich beliebiger Mechnismus möglich bei Anpassung der Ausgansgkonzentrationen und vorhandener Stoffe für Auswertung
#Programm ist für Oberflächenreaktionen 
reaction_mechanism1 ='methane_pox_on_pt.cti'
#reaction_mechanism2 = 'h2o2.cti'
#reaction_mechanism1 ='__MY_methane_pox_on_pt_surface.cti'
               

'''
#####################################################################
# Method 2: Chain of Reactors
#Beipiel der Implemeterung eines pfr_Reaktors wie in den Beispielen pfr.py und surf-pfr.py (hier auch mit surface chemistry)
#Beliebig viele Reaktorschritte, aber entsprechend keine Implementierung von Wänden und damit kein Wärmeaustausch möglich.
#Katalysator befindet sich auf jeder Oberfläche --> Überall kontinierliche Reaktion
#Beipiel dient nur als Referenz, Code für das andere Modell nicht notwendig. 

#####################################################################
# The plug flow reactor is represented by a linear chain of zero-dimensional
# reactors. The gas at the inlet to the first one has the specified inlet
# composition, and for all others the inlet composition is fixed at the
# composition of the reactor immediately upstream. Since in a PFR model there
# is no diffusion, the upstream reactors are not affected by any downstream
# reactors, and therefore the problem may be solved by simply marching from
# the first to last reactor, integrating each one to steady state.
# (This approach is anologous to the one presented in 'surf_pfr.py', which
# additionally includes surface chemistry)

# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.
n_steps = 12

# import the gas model and set the initial conditions
gas2a = ct.Solution(reaction_mechanism1,'gas')

gas2a.TPX = T_0, pressure, composition_0
mass_flow_rate2a = u_0 * gas2a.density * area
dz = length / n_steps
r_vol = area * dz

# import the surface model                                  #surface    
surf = ct.Interface(reaction_mechanism1, 'Pt_surf', [gas2a])   
surf.TP = T_0, ct.one_atm

# catalyst area in one reactor                              #surface
cat_area = cat_area_per_vol * r_vol

# create a new reactor
r2a = ct.IdealGasReactor(gas2a)
r2a.volume = r_vol


# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream2a = ct.Reservoir(gas2a, name='upstream2a')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream2a = ct.Reservoir(gas2a, name='downstream2a')

# Add the reacting surface to the reactor. The area is set to the desired
# catalyst area in the reactor.
rsurf = ct.ReactorSurface(surf, r2a, A=cat_area)                                  

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.
m2a = ct.MassFlowController(upstream2a, r2a, mdot=mass_flow_rate2a)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v2a = ct.PressureController(r2a, downstream2a, master=m2a,K=1e-5) #K=1e-5

sim2a = ct.ReactorNet([r2a])

# define time, space, and other information vectors
z2a = (np.arange(n_steps) + 1) * dz
t_r2a = np.zeros_like(z2a)  # residence time in each reactor
u2a = np.zeros_like(z2a)
t2a = np.zeros_like(z2a)
r2a_density=np.zeros_like(z2a)
states2a = ct.SolutionArray(r2a.thermo)

# iterate through the PFR cells
for n in range(n_steps):
    # Set the state of the reservoir to match that of the previous reactor
    gas2a.TDY = r2a.thermo.TDY
    upstream2a.syncState()
    # integrate the reactor forward in time until steady state is reached
    sim2a.reinitialize()
    sim2a.advance_to_steady_state()
    # compute velocity and transform into time
    u2a[n] = mass_flow_rate2a / area / r2a.thermo.density
    t_r2a[n] = r2a.mass / mass_flow_rate2a  # residence time in this reactor
    t2a[n] = np.sum(t_r2a)
    r2a_density[n]=r2a.thermo.density
    # write output data
    states2a.append(r2a.thermo.state)



plt.figure(1)
plt.plot(z2a, states2a.T,'^-',label='Reactor Chain r2')

#plt.plot(z2, states1_reservoir.T,':', label='Reservoir')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
#plt.show()
plt.savefig('pfr_T_z.png')


plt.figure(2)
plt.plot(t2a, states2a.X[:, gas2a.species_index('CH4')], '^-',label='Reactor Chain R2a')
plt.xlabel('$t$ [s]')
plt.ylabel('$X_{CH_4}$ [-]')
plt.legend(loc=0)


#plt.savefig('pfr_XH2_t.png')

plt.figure(3)
plt.clf()
plt.subplot(3, 1, 1)
h=plt.plot(z2a, states2a.T,'^-')
plt.xlabel('$z$ [m]')
plt.ylabel('Temperature (K)')

plt.subplot(3, 1, 2)
h=plt.plot(z2a, states2a.P,'*-')
plt.xlabel('$z$ [m]')
plt.ylabel('Pressure (Pa)')

plt.subplot(3, 1, 3)
h=plt.plot(z2a, r2a_density,'*-')
plt.xlabel('$z$ [m]')
plt.ylabel('Density(kg/m3)')
#plt.show()
''' 
####END OF METHOD 2 ########################

 
#####################################################################

# Reactor Chain --Single steps 

#####################################################################

# Resolution: The PFR will be simulated by 'n_steps' time steps or by a chain
# of 'n_steps' stirred reactors.
n_steps = 12                          #Hauptprogramm funktioniert nur für diese eine Zahl von Reaktoren. 
#Noch keine Automatisierung implementiert, nur manuelle Anpassung möglich

# import the gas model and set the initial conditions
gas = ct.Solution(reaction_mechanism1, 'gas')
gas.TPX = T_0, pressure, composition_0

mass_flow_rate = u_0 * gas.density * area
dz = length / n_steps
r_vol = area * dz

# import the surface model
# Verwendung einer einzigen Katalysatorfläche scheint nicht zu funktionieren,
# da Temperaturen der Oerfläche und des Reaktors dann nicht übereinstimmen
# surf = ct.Interface(reaction_mechanism1, 'Pt_surf', [gas])
# surf.TPX = T_0, ct.one_atm, initial_coverage_catalyst
#hier haben Reaktoren 5 bis 8 eine Katalysatoroberfläche
surf5 = ct.Interface(reaction_mechanism1, 'Pt_surf', [gas])        
surf6 = ct.Interface(reaction_mechanism1, 'Pt_surf', [gas])        
surf7 = ct.Interface(reaction_mechanism1, 'Pt_surf', [gas])
surf8 = ct.Interface(reaction_mechanism1, 'Pt_surf', [gas]) 
surf5.TP = T_0, ct.one_atm     
surf6.TP = T_0, ct.one_atm
surf7.TP = T_0, ct.one_atm
surf8.TP = T_0, ct.one_atm

# catalyst area in one reactor
cat_area = cat_area_per_vol * r_vol

#Verwendung eines Einzelnen Gases ("gas") statt  ("gas1","gas2" ...) scheint zu funktionieren. ggf. nochmal überprüfen/durchdenken
# create a new reactor
r1 = ct.IdealGasReactor(gas)
r2 = ct.IdealGasReactor(gas)
r3 = ct.IdealGasReactor(gas)
r4 = ct.IdealGasReactor(gas)
r5 = ct.IdealGasReactor(gas)
r6 = ct.IdealGasReactor(gas)
r7 = ct.IdealGasReactor(gas)
r8 = ct.IdealGasReactor(gas)
r9 = ct.IdealGasReactor(gas)
r10 = ct.IdealGasReactor(gas)
r11 = ct.IdealGasReactor(gas)
r12 = ct.IdealGasReactor(gas)

r1.volume = r_vol
r2.volume = r_vol
r3.volume = r_vol
r4.volume = r_vol
r5.volume = r_vol
r6.volume = r_vol
r7.volume = r_vol
r8.volume = r_vol
r9.volume = r_vol
r10.volume = r_vol
r11.volume = r_vol
r12.volume = r_vol


#Wände zwischen Reaktoren für Wärmesaustausch
#Reaktoren mit Katalysator haben hier keine Wand für Wärmeübergang 
if walls==1:
    ct.Wall(r1, r12, A=0.001, U=1.2)    #Wand zwischen Reaktor r1 und r12
    ct.Wall(r2, r11, A=0.001, U=1.2)
    ct.Wall(r3, r10, A=0.001, U=1.2)
    ct.Wall(r4, r9, A=0.001, U=1.2)
    #ct.Wall(r5, r8, A=0.001, U=1.2)
    #ct.Wall(r6, r7, A=0.001, U=1.2)
    print('Walls are active!')

# create a reservoir to represent the reactor immediately upstream. Note
# that the gas object is set already to the state of the upstream reactor
upstream1  = ct.Reservoir(gas, name='upstream1')
upstream2  = ct.Reservoir(gas, name='upstream2')
upstream3  = ct.Reservoir(gas, name='upstream3')
upstream4  = ct.Reservoir(gas, name='upstream4')
upstream5  = ct.Reservoir(gas, name='upstream5')
upstream6  = ct.Reservoir(gas, name='upstream6')
upstream7  = ct.Reservoir(gas, name='upstream7')
upstream8  = ct.Reservoir(gas, name='upstream8')
upstream9  = ct.Reservoir(gas, name='upstream9')
upstream10 = ct.Reservoir(gas, name='upstream10')
upstream11 = ct.Reservoir(gas, name='upstream11')
upstream12 = ct.Reservoir(gas, name='upstream12')

# create a reservoir for the reactor to exhaust into. The composition of
# this reservoir is irrelevant.
downstream1 = ct.Reservoir(gas, name='downstream1')
downstream2 = ct.Reservoir(gas, name='downstream2')
downstream3 = ct.Reservoir(gas, name='downstream3')
downstream4 = ct.Reservoir(gas, name='downstream4')
downstream5 = ct.Reservoir(gas, name='downstream5')
downstream6 = ct.Reservoir(gas, name='downstream6')
downstream7 = ct.Reservoir(gas, name='downstream7')
downstream8 = ct.Reservoir(gas, name='downstream8')
downstream9 = ct.Reservoir(gas, name='downstream9')
downstream10 = ct.Reservoir(gas, name='downstream10')
downstream11 = ct.Reservoir(gas, name='downstream11')
downstream12 = ct.Reservoir(gas, name='downstream12')

# Add the reacting surface to the reactor. The area is set to the desired
# catalyst area in the reactor.

#rsurf1  = ct.ReactorSurface(surf1, r1, A=cat_area)  
#rsurf2  = ct.ReactorSurface(surf2, r2, A=cat_area)  
#rsurf3  = ct.ReactorSurface(surf3, r3, A=cat_area) 
#rsurf4  = ct.ReactorSurface(surf4, r4, A=cat_area) 
rsurf5  = ct.ReactorSurface(surf5, r5, A=cat_area) 
rsurf6  = ct.ReactorSurface(surf6, r6, A=cat_area) 
rsurf7  = ct.ReactorSurface(surf7, r7, A=cat_area) 
rsurf8  = ct.ReactorSurface(surf8, r8, A=cat_area)  
#rsurf9  = ct.ReactorSurface(sur9, r9, A=cat_area) 
#rsurf10 = ct.ReactorSurface(surf10, r10, A=cat_area) 
#rsurf11 = ct.ReactorSurface(surf11, r11, A=cat_area) 
#rsurf12 = ct.ReactorSurface(surf12, r12, A=cat_area)

# The mass flow rate into the reactor will be fixed by using a
# MassFlowController object.

m1 = ct.MassFlowController(upstream1, r1, mdot=mass_flow_rate)
m2 = ct.MassFlowController(upstream2, r2, mdot=mass_flow_rate)
m3 = ct.MassFlowController(upstream3, r3, mdot=mass_flow_rate)
m4 = ct.MassFlowController(upstream4, r4, mdot=mass_flow_rate)
m5 = ct.MassFlowController(upstream5, r5, mdot=mass_flow_rate)
m6 = ct.MassFlowController(upstream6, r6, mdot=mass_flow_rate)
m7 = ct.MassFlowController(upstream7, r7, mdot=mass_flow_rate)
m8 = ct.MassFlowController(upstream8, r8, mdot=mass_flow_rate)
m9 = ct.MassFlowController(upstream9, r9, mdot=mass_flow_rate)
m10 = ct.MassFlowController(upstream10, r10, mdot=mass_flow_rate)
m11 = ct.MassFlowController(upstream11, r11, mdot=mass_flow_rate)
m12 = ct.MassFlowController(upstream12, r12, mdot=mass_flow_rate)

# We need an outlet to the downstream reservoir. This will determine the
# pressure in the reactor. The value of K will only affect the transient
# pressure difference.
v1 = ct.PressureController(r1, downstream1, master=m1, K=1e-5)
v2 = ct.PressureController(r2, downstream2, master=m2, K=1e-5)
v3 = ct.PressureController(r3, downstream3, master=m3, K=1e-5)
v4 = ct.PressureController(r4, downstream4, master=m4, K=1e-5)
v5 = ct.PressureController(r5, downstream5, master=m5, K=1e-5)
v6 = ct.PressureController(r6, downstream6, master=m6, K=1e-5)
v7 = ct.PressureController(r7, downstream7, master=m7, K=1e-5)
v8 = ct.PressureController(r8, downstream8, master=m8, K=1e-5)
v9 = ct.PressureController(r9, downstream9, master=m9, K=1e-5)
v10 = ct.PressureController(r10, downstream10, master=m10, K=1e-5)
v11 = ct.PressureController(r11, downstream11, master=m11, K=1e-5)
v12 = ct.PressureController(r12, downstream12, master=m12, K=1e-5)


sim2 = ct.ReactorNet([r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12]) #alle Reaktoren werden in jedem Schritt berechnet

# Define time, space, and other information vectors

# "Länge" des Gesamtreaktors nach jedem Teilreaktor"
z_vec=[]
for i in range(1,n_steps+1):
    z_vec.append(i*dz)

#Zustände in den Reaktoren
states1 = ct.SolutionArray(r1.thermo)
states2 = ct.SolutionArray(r2.thermo)
states3 = ct.SolutionArray(r3.thermo)
states4 = ct.SolutionArray(r4.thermo)
states5 = ct.SolutionArray(r5.thermo)
states6 = ct.SolutionArray(r6.thermo)
states7 = ct.SolutionArray(r7.thermo)
states8 = ct.SolutionArray(r8.thermo)
states9 = ct.SolutionArray(r9.thermo)
states10 = ct.SolutionArray(r10.thermo)
states11 = ct.SolutionArray(r11.thermo)
states12 = ct.SolutionArray(r12.thermo)

#Exttra Temperaturvektor zur einfacheren Auswertung der Temperaturen
r1_T_vec=[]
r2_T_vec=[]
r3_T_vec=[]
r4_T_vec=[]
r5_T_vec=[]
r6_T_vec=[]
r7_T_vec=[]
r8_T_vec=[]
r9_T_vec=[]
r10_T_vec=[]
r11_T_vec=[]
r12_T_vec=[]

#Vektoren zur Auswertung
temp_profile=[]
p_profile=[]
X_CH4_profile=[]
X_O2_profile=[]
X_CO2_profile=[]
X_H2O_profile=[]

# iterate through the PFR cells

for n in range(steps):

    #Reaktor 1
    
    # Gas nach jededem Durchlauf auf Ausgangszusammensetzung zurücksetzen  
    gas.TPX = T_0, pressure, composition_0

    # Set the state of the reservoir to match that of the previous reactor
    upstream1.syncState()                   
 
    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()        

    # compute velocity and transform into time
    u1 = mass_flow_rate / area / r1.thermo.density
    t_r1 = r1.mass / mass_flow_rate  # residence time in this reactor
    t1 = t_r1
    r1_density=r1.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)

    # Übergang in Reaktor 2 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r1.thermo.TDY
    upstream2.syncState() 
    
    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    u2= mass_flow_rate / area / r2.thermo.density
    t_r2 = r2.mass / mass_flow_rate  # residence time in this reactor
    t2 = t1+t_r2
    r2_density=r2.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    # Übergang in Reaktor 3 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r2.thermo.TDY 
    upstream3.syncState()   

    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()          #step 3 Durchgnag 3
   
    # compute velocity and transform into time
    u3 = mass_flow_rate / area / r3.thermo.density
    t_r3 = r3.mass / mass_flow_rate  # residence time in this reactor
    t3 = t2+t_r3
    r3_density=r3.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    # Übergang in Reaktor 4 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r3.thermo.TDY
    upstream4.syncState() 
    
    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    u4 = mass_flow_rate / area / r4.thermo.density
    t_r4 = r4.mass / mass_flow_rate  # residence time in this reactor
    t4 = t3+t_r4
    r4_density=r4.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    # Übergang in Reaktor 5 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r4.thermo.TDY
    upstream5.syncState()

    #Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    # surf5_tmp=surf5.TPX
    #print("surf5_tmp =",surf5_tmp[0],r5.T)

    #integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    #Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    #surf5_tmp=surf5.TPX
    #print("surf5_tmp =",surf5_tmp[0],r5.T)

    # compute velocity and transform into time
    u5 = mass_flow_rate / area / r5.thermo.density
    t_r5 = r5.mass / mass_flow_rate  # residence time in this reactor
    t5 = t4+t_r5
    r5_density=r5.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)

    #####################################################################
    # Übergang in Reaktor 6 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r5.thermo.TDY
    upstream6.syncState()

    # Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    # surf6_tmp=surf6.TPX
    # print("surf6_tmp =",surf6_tmp[0],r6.T)

    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    # Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    # surf6_tmp=surf6.TPX
    # print("surf6_tmp =",surf6_tmp[0],r6.T)

    # compute velocity and transform into time
    u6 = mass_flow_rate / area / r6.thermo.density
    t_r6 = r6.mass / mass_flow_rate  # residence time in this reactor
    t6 = t5+t_r6
    r6_density=r6.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Übergang in Reaktor 7 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r6.thermo.TDY
    upstream7.syncState()

    # Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    # surf7_tmp=surf7.TPX
    # print("surf7_tmp =",surf7_tmp[0],r7.T)

    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    # Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    # surf7_tmp=surf7.TPX
    # print("surf7_tmp =",surf7_tmp[0],r7.T)

    # compute velocity and transform into time
    u7 = mass_flow_rate / area / r7.thermo.density
    t_r7 = r7.mass / mass_flow_rate  # residence time in this reactor
    t7 = t6+t_r7
    r7_density=r7.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)

    # Übergang in Reaktor 8 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r7.thermo.TDY
    upstream8.syncState()

    # Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    # surf8_tmp=surf8.TPX
    # print("surf8_tmp =",surf8_tmp[0],r8.T)

    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    # Kontrolle der Temperaturen in Reaktor und Oberfläche (provisorisch)
    # surf8_tmp=surf8.TPX
    # print("surf8_tmp =",surf8_tmp[0],r8.T)

    # compute velocity and transform into time
    u8 = mass_flow_rate / area / r8.thermo.density
    t_r8 = r8.mass / mass_flow_rate  # residence time in this reactor
    t8 = t7+t_r8
    r8_density=r8.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)

    # Übergang in Reaktor 9 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r8.thermo.TDY
    upstream9.syncState()

    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    # compute velocity and transform into time
    u9 = mass_flow_rate / area / r9.thermo.density
    t_r9 = r9.mass / mass_flow_rate  # residence time in this reactor
    t9 = t8+t_r9
    r9_density=r9.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)

    # Übergang in Reaktor 10 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r9.thermo.TDY
    upstream10.syncState()

    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    # compute velocity and transform into time
    u10 = mass_flow_rate / area / r6.thermo.density
    t_r10 = r10.mass / mass_flow_rate  # residence time in this reactor
    t10 = t9+t_r10
    r10_density=r10.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)

    # Übergang in Reaktor 11 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r10.thermo.TDY
    upstream11.syncState()

    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()

    # compute velocity and transform into time
    u11 = mass_flow_rate / area / r11.thermo.density
    t_r11 = r11.mass / mass_flow_rate  # residence time in this reactor
    t11 = t10+t_r11
    r11_density=r11.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)


    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    # Übergang in Reaktor 12 &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # Set the state of the reservoir to match that of the previous reactor
    gas.TDY = r11.thermo.TDY
    upstream12.syncState()
    
    # integrate the reactor forward in time until steady state is reached
    sim2.reinitialize()                     
    sim2.advance_to_steady_state()          #step 3 Durchgnag 3

    # compute velocity and transform into time
    u12 = mass_flow_rate / area / r12.thermo.density
    t_r12 = r12.mass / mass_flow_rate  # residence time in this reactor
    t12 = t11+t_r12
    r12_density=r12.thermo.density

    #Abspeichern der Temperaturen in jedem Reaktor
    #Abspeichern der Temperaturen in jedem Reaktor
    r1_T_vec.append(r1.T)
    r2_T_vec.append(r2.T)
    r3_T_vec.append(r3.T)
    r4_T_vec.append(r4.T)
    r5_T_vec.append(r5.T)
    r6_T_vec.append(r6.T)
    r7_T_vec.append(r7.T)
    r8_T_vec.append(r8.T)
    r9_T_vec.append(r9.T)
    r10_T_vec.append(r10.T)
    r11_T_vec.append(r11.T)
    r12_T_vec.append(r12.T)

    #Abspeichern aller thermo.state Variablen in jedem Reaktor
    states1.append(r1.thermo.state)
    states2.append(r2.thermo.state)
    states3.append(r3.thermo.state)
    states4.append(r4.thermo.state)
    states5.append(r5.thermo.state)
    states6.append(r6.thermo.state)
    states7.append(r7.thermo.state)
    states8.append(r8.thermo.state)
    states9.append(r9.thermo.state)
    states10.append(r10.thermo.state)
    states11.append(r11.thermo.state)
    states12.append(r12.thermo.state)
    
    #Abspeichern des Temperaturprofils über alle Stützpunkte nach jedem vollen Durchlauf 
    temp_profile.append(r1.T) 
    temp_profile.append(r2.T)
    temp_profile.append(r3.T) 
    temp_profile.append(r4.T) 
    temp_profile.append(r5.T) 
    temp_profile.append(r6.T) 
    temp_profile.append(r7.T) 
    temp_profile.append(r8.T)
    temp_profile.append(r9.T)
    temp_profile.append(r10.T)
    temp_profile.append(r11.T)
    temp_profile.append(r12.T)

    #Abspeichern des Druckprofils über alle Stützpunkte nach jedem vollen Durchlauf 
    p_profile.append(r1.thermo.P)
    p_profile.append(r2.thermo.P)
    p_profile.append(r3.thermo.P)
    p_profile.append(r4.thermo.P)
    p_profile.append(r5.thermo.P)
    p_profile.append(r6.thermo.P)
    p_profile.append(r7.thermo.P)
    p_profile.append(r8.thermo.P)
    p_profile.append(r9.thermo.P)
    p_profile.append(r10.thermo.P)
    p_profile.append(r11.thermo.P)
    p_profile.append(r12.thermo.P)

    #Abspeichern einiger Konzentrationsprofile über alle Stützpunkte nach jedem vollen Durchlauf 
    X_CH4_profile.append(r1.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r2.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r3.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r4.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r5.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r6.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r7.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r8.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r9.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r10.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r11.thermo.X[gas.species_index('CH4')])
    X_CH4_profile.append(r12.thermo.X[gas.species_index('CH4')])

    
    X_O2_profile.append(r1.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r2.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r3.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r4.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r5.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r6.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r7.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r8.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r9.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r10.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r11.thermo.X[gas.species_index('O2')])
    X_O2_profile.append(r12.thermo.X[gas.species_index('O2')])

    X_CO2_profile.append(r1.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r2.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r3.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r4.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r5.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r6.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r7.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r8.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r9.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r10.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r11.thermo.X[gas.species_index('CO2')])
    X_CO2_profile.append(r12.thermo.X[gas.species_index('CO2')])

    X_H2O_profile.append(r1.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r2.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r3.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r4.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r5.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r6.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r7.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r8.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r9.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r10.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r11.thermo.X[gas.species_index('H2O')])
    X_H2O_profile.append(r12.thermo.X[gas.species_index('H2O')])

    #Ausgabe der Reaktordurchläufe zur Kontrolle
    print('Reaktordurchlauf:\t', n+1)

#Ausagabe der Zusammenstzung am Austritt
gas.TDY = r12.thermo.TDY
gas()

#####################################################################
# Compare Results in matplotlib
#####################################################################

#Plots 4 bis 9 zeigen beispielhaft Temperaturentwicklung innerhalb eines Reaktors bei Auswertung für jeden Berechnungsschritt
#hilfreich um Prinzip der Schritweisen weitergabe des Reaktionsgemischs und anschließende Wärmeübertragung durch die Wand zu verstehen
#ansonsten eher unwichtig
#Länge z in der x-Achse könnte auch durch  Reaktoren 1 bis 12 ersetzt werden und stehte für die aktuelle Position des Reaktionsgemischs

reactor_steps=np.arange(1,12+1)

plt.figure(4)
for n in range(steps):
    plt.plot(z_vec, r1_T_vec[0+(n_steps)*n:(n_steps)+(n_steps)*n],'^-',label='T of r1')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]') 
plt.legend(loc=0)
plt.title('Reactor 1')

'''
plt.figure(5)
for n in range(steps):
    plt.plot(z_vec, r2_T_vec[0+(n_steps)*n:(n_steps)+(n_steps)*n],'^-',label='T of r2')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.title('Reactor 2')

plt.figure(6)
for n in range(steps):
    plt.plot(z_vec, r3_T_vec[0+n_steps*n:n_steps+n_steps*n],'^-',label='T of r3')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.title('Reactor 3')

plt.figure(7)
for n in range(steps):
    plt.plot(z_vec, r4_T_vec[0+n_steps*n:n_steps+n_steps*n],'^-',label='T of r4')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.title('Reactor 4')

plt.figure(8)
for n in range(steps):
    plt.plot(z_vec, r5_T_vec[0+n_steps*n:n_steps+n_steps*n],'^-',label='T of r5')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.title('Reactor 5')

plt.figure(9)
for n in range(steps):
    plt.plot(z_vec, r6_T_vec[0+n_steps*n:n_steps+n_steps*n],'^-',label='T of r6')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.title('Reactor 6')

plt.figure(10)
for n in range(steps):
    plt.plot(z_vec, r7_T_vec[0+n_steps*n:n_steps+n_steps*n],'^-',label='T of r7')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.title('Reactor 7')

plt.figure(11)
for n in range(steps):
    plt.plot(z_vec, r12_T_vec[0+n_steps*n:n_steps+n_steps*n],'^-',label='T of r12')
plt.xlabel('$z$ [m]')
plt.ylabel('$T$ [K]')
plt.legend(loc=0)
plt.title('Reactor 12')
'''
# Plots: 12-16 "full timeline". Alle Werte werden nach jedem Teilschritt (Weitergabe der Gase zum nächsten Reaktor) in allen 12 Reaktoren abgespeichert und hier aufgetragen
# Nach jeweils 12 Schritten erfolget ein neuer Durchlauf der 12 Reaktoren mit frischen Edukten 
# x Achsen anzahl der Gsamtschritte also Durchläufe*12
plt.figure(12)
plt.plot(np.arange(2,len(states1.T)+2),states1.T,'^-',label='T of r1')
plt.plot(np.arange(2,len(states2.T)+2),states2.T,'^-',label='T of r2')
plt.plot(np.arange(2,len(states3.T)+2),states3.T,'^-',label='T of r3')
plt.plot(np.arange(2,len(states4.T)+2),states4.T,'^-',label='T of r4')
plt.plot(np.arange(2,len(states5.T)+2),states5.T,'^-',label='T of r5')
plt.plot(np.arange(2,len(states6.T)+2),states6.T,'^-',label='T of r6')
plt.plot(np.arange(2,len(states7.T)+2),states7.T,'^-',label='T of r7')
plt.plot(np.arange(2,len(states8.T)+2),states8.T,'^-',label='T of r8')
plt.plot(np.arange(2,len(states9.T)+2),states9.T,'^-',label='T of r9')
plt.plot(np.arange(2,len(states10.T)+2),states10.T,'^-',label='T of r10')
plt.plot(np.arange(2,len(states11.T)+2),states11.T,'^-',label='T of r11')
plt.plot(np.arange(2,len(states12.T)+2),states12.T,'^-',label='T of r12')
plt.xlabel('$Reaktorschritte$')
plt.ylabel('$T$ [K]')
plt.title('Temperaturverläufe in alle Reaktoren')
plt.legend(loc=0)

plt.figure(13)
plt.plot(np.arange(2,len(states2.X[:, gas.species_index('CH4')])+2),states2.X[:, gas.species_index('CH4')],'^-',label='CH4')
plt.plot(np.arange(2,len(states2.X[:, gas.species_index('O2')])+2),states2.X[:, gas.species_index('O2')],'^-',label='O2')
plt.plot(np.arange(2,len(states2.X[:, gas.species_index('CO2')])+2),states2.X[:, gas.species_index('CO2')],'^-',label='CO2')
plt.plot(np.arange(2,len(states2.X[:, gas.species_index('H2O')])+2),states2.X[:, gas.species_index('H2O')],'^-',label='H2O')
plt.legend(loc=0)
plt.xlabel('$Reaktorschritte$')
plt.ylabel('$X$ [-]') #Molefraction
plt.title('Molefrac in Reactor 2')

plt.figure(14)
plt.plot(np.arange(2,len(states5.X[:, gas.species_index('CH4')])+2),states5.X[:, gas.species_index('CH4')],'^-',label='CH4')
plt.plot(np.arange(2,len(states5.X[:, gas.species_index('O2')])+2),states5.X[:, gas.species_index('O2')],'^-',label='O2')
plt.plot(np.arange(2,len(states5.X[:, gas.species_index('CO2')])+2),states5.X[:, gas.species_index('CO2')],'^-',label='CO2')
plt.plot(np.arange(2,len(states5.X[:, gas.species_index('H2O')])+2),states5.X[:, gas.species_index('H2O')],'^-',label='H2O')
plt.legend(loc=0)
plt.xlabel('$Reaktorschritte$')
plt.ylabel('$X$ [-]') #Molefraction
plt.title('Molefrac in Reactor 5')

'''
plt.figure(15)
plt.plot(np.arange(2,len(states6.X[:, gas.species_index('CH4')])+2),states6.X[:, gas.species_index('CH4')],'^-',label='CH4')
plt.plot(np.arange(2,len(states6.X[:, gas.species_index('O2')])+2),states6.X[:, gas.species_index('O2')],'^-',label='O2')
plt.plot(np.arange(2,len(states6.X[:, gas.species_index('CO2')])+2),states6.X[:, gas.species_index('CO2')],'^-',label='CO2')
plt.plot(np.arange(2,len(states6.X[:, gas.species_index('H2O')])+2),states6.X[:, gas.species_index('H2O')],'^-',label='H2O')
plt.legend(loc=0)
plt.xlabel('$Reaktorschritte$')
plt.ylabel('$X$ [-]') #Molefraction
plt.title('Molefrac in Reactor 6')
'''

plt.figure(16)
plt.plot(np.arange(2,len(states12.X[:, gas.species_index('CH4')])+2),states12.X[:, gas.species_index('CH4')],'^-',label='CH4')
plt.plot(np.arange(2,len(states12.X[:, gas.species_index('O2')])+2),states12.X[:, gas.species_index('O2')],'^-',label='O2')
plt.plot(np.arange(2,len(states12.X[:, gas.species_index('CO2')])+2),states12.X[:, gas.species_index('CO2')],'^-',label='CO2')
plt.plot(np.arange(2,len(states12.X[:, gas.species_index('H2O')])+2),states12.X[:, gas.species_index('H2O')],'^-',label='H2O')
plt.legend(loc=0)
plt.xlabel('$Reaktorschritte$')
plt.ylabel('$X$ [-]') #Molefraction
plt.title('Molefrac in Reactor 12')


# Profil im Reaktor nach jedem vollen Durchlauf. Werte in allen 12 Reaktoren werden jeweils nach der Reaktion im letzten Reaktor (Reaktor 12) abgespeichert.
# Nach jedem vollen Durchlauf wird das Pofil mit den Sützpunkten Reaktor 1 bis 12 geplottet.
#  -> Itertiver Verlauf der Profie im Gesamtreaktor bis zu stationärem Bestriebszustand

plt.figure(20)
#Temperatur
reactor_steps=np.arange(1,12+1)

for n in range(steps):
    plt.plot(reactor_steps,temp_profile[0+n_steps*n:n_steps+n_steps*n],'^-',label='Temp Profile')
plt.xlabel('$Reaktor$')
plt.ylabel('$T$ [K]')
plt.title('Entwicklung des Temperaturprofils')
plt.legend(loc=0)

plt.figure(21)
for n in range(steps):
    plt.plot(reactor_steps,p_profile[0+n_steps*n:n_steps+n_steps*n],'^-',label='Pressure Profile')
plt.xlabel('$Reaktor$')
plt.ylabel('$P$ [Pa]') 
plt.title('Entwicklung des Druckprofils')
plt.legend(loc=0)

plt.figure(22)
for n in range(steps):
    plt.plot(reactor_steps,X_CH4_profile[0+n_steps*n:n_steps+n_steps*n],'^-',label='X CH4 Profile')
plt.xlabel('$Reaktor$')
plt.ylabel('$X$ [-]') 
plt.title('Entwicklung des CH4 Profils')
plt.legend(loc=0)

plt.figure(23)
for n in range(steps):
    plt.plot(reactor_steps,X_O2_profile[0+n_steps*n:n_steps+n_steps*n],'^-',label='X O2 Profile')
plt.xlabel('$Reaktor$')
plt.ylabel('$X$ [-]') 
plt.title('Entwicklung des O2 Profils')
plt.legend(loc=0)


plt.figure(24)
for n in range(steps):
    plt.plot(reactor_steps,X_CO2_profile[0+n_steps*n:n_steps+n_steps*n],'^-',label='X CO2 Profile')
plt.xlabel('$Reaktor$')
plt.ylabel('$X$ [-]') 
plt.title('Entwicklung des CO2 Profils')
plt.legend(loc=0)

plt.figure(25)
for n in range(steps):
    plt.plot(reactor_steps,X_H2O_profile[0+n_steps*n:n_steps+n_steps*n],'^-',label='X H2O Profile')
plt.xlabel('$Reaktor$')
plt.ylabel('$X$ [-]') 
plt.title('Entwicklung des H20 Profils')
plt.legend(loc=0)

plt.show()