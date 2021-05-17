import cantera as ct
import numpy as np

gas = ct.Solution('mech_13.yaml')
gas.TPX = 870. , ct.one_atm , 'C3H8:10, H2:1'
gas()