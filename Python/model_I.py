# Zooplankton community model I
# Two unstructured consumer species competing for two resources
# For units and references, see Appendix S2, Table 2
# Created by Wojciech Uszko (2020)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Body masses (ng dry weight):

B_C1 = 100     # consumer 1
B_C2 = 1000    # consumer 2

# Temperature- or body mass-independent parameters:

deltaRS = 0.1  # resource 1 supply rate
deltaRL = 0.1  # resource 2 supply rate

q       = 0    # functional response (Hill) exponent; if =0 then type II

p       = 0.85 # diet preference
pC1RS   = p
pC1RL   = 1-pC1RS
pC2RS   = 1-pC1RS
pC2RL   = pC1RS

betaC1  = 0.6  # consumer 1 conversion efficiency
betaC2  = 0.6  # consumer 2 conversion efficiency

HC1RS   = 0.2  # half-saturation constant
HC1RL   = 0.2  # half-saturation constant
HC2RS   = 0.2  # half-saturation constant
HC2RL   = 0.2  # half-saturation constant

muC1    = 0.01 # consumer 1 background mortality rate
muC2    = 0.01 # consumer 2 background mortality rate

# Ambient temperature (Kelvin):

T = 273.15 + 20

"""
# Temperature- or body mass-dependent parameters
# Without size-temperature interaction:

# Resource supply density:
RSmax    = 0.0042 * np.exp( 0.151/(0.00008617*T) )
RLmax    = RSmax
 
# Consumer maximum ingestion rate:
IC1RSmax = (19 * (B_C1**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C1
IC1RLmax = (19 * (B_C1**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C1
IC2RSmax = (19 * (B_C2**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C2
IC2RLmax = (19 * (B_C2**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C2

# Consumer metabolic rate:
mC1      = (850000000 * (B_C1**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_C1
mC2      = (850000000 * (B_C2**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_C2
"""


# Temperature- or body mass-dependent parameters
# With size-temperature interaction in Rmax and in temperature optimum of Imax:

# Resource supply density:
RSmax    = 0.0042 * np.exp( 0.151/(0.00008617*T) )
RLmax    = (5.88* 10**(-7)) * np.exp( 0.37564/(0.00008617*T) )

# Consumer maximum ingestion rate:
IC1RSmax = (19 * (B_C1**(0.7)) * np.exp(-((T-(273.15+24))**2)/(2*(8**2)))) / B_C1
IC1RLmax = (19 * (B_C1**(0.7)) * np.exp(-((T-(273.15+24))**2)/(2*(8**2)))) / B_C1
IC2RSmax = (19 * (B_C2**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C2
IC2RLmax = (19 * (B_C2**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C2

# Consumer metabolic rate:
mC1      = (850000000 * (B_C1**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_C1
mC2      = (850000000 * (B_C2**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_C2


"""
# Temperature- or body mass-dependent parameters
# With size-temperature interaction in Rmax and in metabolic rate:

# Resource supply density:
RSmax    = 0.0042 * np.exp( 0.151/(0.00008617*T) )
RLmax    = (5.88* 10**(-7)) * np.exp( 0.37564/(0.00008617*T) )

# Consumer maximum ingestion rate:
IC1RSmax = (19 * (B_C1**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C1
IC1RLmax = (19 * (B_C1**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C1
IC2RSmax = (19 * (B_C2**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C2
IC2RLmax = (19 * (B_C2**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_C2

# Consumer metabolic rate:
mC1      = (850000000 * (B_C1**(0.7 + 0.0005*T)) * np.exp( -0.56/(0.00008617*T) )) / B_C1
mC2      = (850000000 * (B_C2**(0.7 + 0.0005*T)) * np.exp( -0.56/(0.00008617*T) )) / B_C2
"""

# Specify the model:
def model(X,t):
    # Variables:
    RS = X[0]   # small resource biomass density
    RL = X[1]   # large resource biomass density
    C1 = X[2]   # consumer 1 biomass density
    C2 = X[3]   # consumer 2 biomass density
    # Ingestion rates:
    IC1RS = ( ( pC1RS * (IC1RSmax/(HC1RS**(1+q))) * RS**(1+q) + 0 * (IC1RLmax/(HC1RL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pC1RS/(HC1RS**(1+q))) * RS**(1+q) + (pC1RL/(HC1RL**(1+q))) * RL**(1+q) ) )
    IC1RL = ( ( 0 * (IC1RSmax/(HC1RS**(1+q))) * RS**(1+q) + pC1RL * (IC1RLmax/(HC1RL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pC1RS/(HC1RS**(1+q))) * RS**(1+q) + (pC1RL/(HC1RL**(1+q))) * RL**(1+q) ) )
    IC2RS = ( ( pC2RS * (IC2RSmax/(HC2RS**(1+q))) * RS**(1+q) + 0 * (IC2RLmax/(HC2RL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pC2RS/(HC2RS**(1+q))) * RS**(1+q) + (pC2RL/(HC2RL**(1+q))) * RL**(1+q) ) )
    IC2RL = ( ( 0 * (IC2RSmax/(HC2RS**(1+q))) * RS**(1+q) + pC2RL * (IC2RLmax/(HC2RL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pC2RS/(HC2RS**(1+q))) * RS**(1+q) + (pC2RL/(HC2RL**(1+q))) * RL**(1+q) ) )
    # ODE system:
    dRSdt = deltaRS*(RSmax - RS) - IC1RS*C1 - IC2RS*C2
    dRLdt = deltaRL*(RLmax - RL) - IC1RL*C1 - IC2RL*C2
    dC1dt =  betaC1*(IC1RS+IC1RL)*C1 - mC1*C1 - muC1*C1
    dC2dt =  betaC2*(IC2RS+IC2RL)*C2 - mC2*C2 - muC2*C2
    return np.array([dRSdt, dRLdt, dC1dt, dC2dt])

# Initial densities for RS, RL, C1, C2
X0 = np.array([0.01, 0.01, 0.01, 0.01])

# Time range
t = np.linspace(0,300,1000)

# Solve ODE
X = odeint(model,X0,t)

# Plot results
RS,RL,C1,C2 = np.transpose(X)

plt.figure()
plt.plot(t, RS, 'g-', label='RS', linewidth=1.0)
plt.plot(t, RL, 'g-', label='RL', linewidth=2.5)
plt.legend(loc='upper right')
plt.xlabel('Time (day)')
plt.ylabel('Density (mg/L)')
plt.show()

plt.figure()
plt.plot(t, C1, 'k-', label='C1', linewidth=1.0)
plt.plot(t, C2, 'k-', label='C2', linewidth=2.5)
plt.legend(loc='upper right')
plt.xlabel('Time (day)')
plt.ylabel('Density (mg/L)')
plt.show()