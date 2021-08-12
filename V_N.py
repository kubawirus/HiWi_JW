import cantera as ct
def vol_norm(v1, gas, T_0, p):

    q1= ct.Quantity(gas, mass = gas.density*v1)
    q1.TP = T_0, p
    v2 = q1.volume
    print(v2)

    return v2


