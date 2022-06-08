from cantera import *
import cantera as ct
import math

# combustion chamber volume and mass flows assessment based on training manual of Rocketdyne F-1 [1]

#lambda_air = 1.5  # air–fuel equivalence ratio
from cantera import Valve

oxy = ct.Solution('gri30.cti')
oxy.TPY = 800, 1.2797e+7, 'O2:1.0'  # according to literature
oxy_in = ct.Reservoir(oxy)
oxy_mdot = ct.Quantity(oxy, mass=1789)

fuel = ct.Solution('Dagaut_Ori.cti')
fuel.TPY = 300, 3e+05, 'NC10H22:0.74,PHC3H7:0.15,CYC9H18:0.11'
fuel_in = ct.Reservoir(fuel)
fuel_mdot = ct.Quantity(fuel, mass=788)
fuel_mdot.mass = oxy_mdot.mass / 2.27

# igniter (like in "combustor.py")
fuel.TPX = 1500, 2e+06, 'H:1.0'
igniter = ct.Reservoir(fuel)

fuel.TPX = 2100, 1.2e+6, 'N2:1.0'  # combustion chamber already hot, otherwise some mechanims don't integrate properly when combustor temperature is below 1000 K
combustor = ct.IdealGasReactor(fuel, energy='on')
combustor.volume = 0.7

# exhaust reservoir
fuel.TPX = 300, ct.one_atm, 'N2:1.0'
exhaust = ct.Reservoir(fuel)
m1 = ct.MassFlowController(fuel_in, combustor, mdot=fuel_mdot.mass)
m2 = ct.MassFlowController(oxy_in, combustor, mdot=oxy_mdot.mass)

fwhm = 0.01
amplitude = 0.1
t0 = 0.2
igniter_mdot = lambda t: amplitude * math.exp(-(t - t0) ** 2 * 4 * math.log(2) / fwhm ** 2)
m3 = ct.MassFlowController(igniter, combustor, mdot=igniter_mdot)

v = ct.Valve(combustor, exhaust, K=0.00037)

sim = ct.ReactorNet([combustor])

time = 0
Tprev = combustor.T

states = ct.SolutionArray(fuel, extra=['t', 'tres'])

for n in range(5000):
    time += 0.0002
    sim.advance(time)
    #tres = combustor.mass / v.mdot(time)
    Tnow = combustor.T
    states.append(fuel.state, t=time, tres=1)

states.write_csv('combustor.csv', cols=('t', 'T', 'P', 'X'))
print(combustor.thermo.report())

#Specific Impulse
c_p=2114.5 #J/(kg*K)
c_v=1727.5 #J/(kg*K)
k=c_p/c_v
R=8314.46/21.488
Isp=((k*R*3691.6/(k-1)*(1-(ct.one_atm/(7.0665*10**6))**((k-1)/k)))**(1/2))/9.81

print(Isp)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(states.t, states.T)
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.title('Temperature')
plt.savefig('T.pdf')
plt.show()

plt.figure()
plt.plot(states.t, states.P)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.tight_layout()
plt.title('Pressure')
plt.savefig('P.pdf')
plt.show()
