# Zooplankton community model II
# One stage-structured consumer species feeding on two resources
# For units and references, see Appendix S2, Table 2
# Created by Wojciech Uszko (2020)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Body masses (ng dry weight):

B_J = 100     # juvenile
B_A = 1000    # adult

# Temperature- or body mass-independent parameters:

deltaRS = 0.1  # resource 1 supply rate
deltaRL = 0.1  # resource 2 supply rate

q       = 0    # functional response (Hill) exponent; if =0 then type II

p       = 0.85 # diet preference
pJRS    = p
pJRL    = 1-pJRS
pARS    = 1-pJRS
pARL    = pJRS

betaJ  = 0.6  # juvenile conversion efficiency
betaA  = 0.6  # adult conversion efficiency

HJRS   = 0.2  # half-saturation constant
HJRL   = 0.2  # half-saturation constant
HARS   = 0.2  # half-saturation constant
HARL   = 0.2  # half-saturation constant

zJA     = 0.1 # juvenile-to-adult mass ratio

muJ    = 0.01 # juvenile background mortality rate
muA    = 0.01 # adult background mortality rate

# Ambient temperature (Kelvin):

T = 273.15 + 20

"""
# Temperature- or body mass-dependent parameters
# Without size-temperature interaction:

# Resource supply density:
RSmax    = 0.0042 * np.exp( 0.151/(0.00008617*T) )
RLmax    = RSmax
 
# Consumer maximum ingestion rate:
IJRSmax = (19 * (B_J**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_J
IJRLmax = (19 * (B_J**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_J
IARSmax = (19 * (B_A**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_A
IARLmax = (19 * (B_A**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_A

# Consumer metabolic rate:
mJ      = (850000000 * (B_J**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_J
mA      = (850000000 * (B_A**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_A
"""


# Temperature- or body mass-dependent parameters
# With size-temperature interaction in Rmax and in temperature optimum of Imax:

# Resource supply density:
RSmax    = 0.0042 * np.exp( 0.151/(0.00008617*T) )
RLmax    = (5.88* 10**(-7)) * np.exp( 0.37564/(0.00008617*T) )

# Consumer maximum ingestion rate:
IJRSmax = (19 * (B_J**(0.7)) * np.exp(-((T-(273.15+24))**2)/(2*(8**2)))) / B_J
IJRLmax = (19 * (B_J**(0.7)) * np.exp(-((T-(273.15+24))**2)/(2*(8**2)))) / B_J
IARSmax = (19 * (B_A**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_A
IARLmax = (19 * (B_A**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_A

# Consumer metabolic rate:
mJ      = (850000000 * (B_J**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_J
mA      = (850000000 * (B_A**0.7) * np.exp( -0.56/(0.00008617*T) )) / B_A


"""
# Temperature- or body mass-dependent parameters
# With size-temperature interaction in Rmax and in metabolic rate:

# Resource supply density:
RSmax    = 0.0042 * np.exp( 0.151/(0.00008617*T) )
RLmax    = (5.88* 10**(-7)) * np.exp( 0.37564/(0.00008617*T) )

# Consumer maximum ingestion rate:
IJRSmax = (19 * (B_J**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_J
IJRLmax = (19 * (B_J**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_J
IARSmax = (19 * (B_A**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_A
IARLmax = (19 * (B_A**(0.7)) * np.exp(-((T-(273.15+20))**2)/(2*(8**2)))) / B_A

# Consumer metabolic rate:
mJ      = (850000000 * (B_J**(0.7 + 0.0005*T)) * np.exp( -0.56/(0.00008617*T) )) / B_J
mA      = (850000000 * (B_A**(0.7 + 0.0005*T)) * np.exp( -0.56/(0.00008617*T) )) / B_A
"""

# Specify the model:
def model(X,t):
    # Variables:
    RS = X[0]   # small resource biomass density
    RL = X[1]   # large resource biomass density
    J = X[2]    # juvenile biomass density
    A = X[3]    # adult biomass density
    # Ingestion rates:
    IJRS = ( ( pJRS * (IJRSmax/(HJRS**(1+q))) * RS**(1+q) + 0 * (IJRLmax/(HJRL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pJRS/(HJRS**(1+q))) * RS**(1+q) + (pJRL/(HJRL**(1+q))) * RL**(1+q) ) )
    IJRL = ( ( 0 * (IJRSmax/(HJRS**(1+q))) * RS**(1+q) + pJRL * (IJRLmax/(HJRL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pJRS/(HJRS**(1+q))) * RS**(1+q) + (pJRL/(HJRL**(1+q))) * RL**(1+q) ) )
    IARS = ( ( pARS * (IARSmax/(HARS**(1+q))) * RS**(1+q) + 0 * (IARLmax/(HARL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pARS/(HARS**(1+q))) * RS**(1+q) + (pARL/(HARL**(1+q))) * RL**(1+q) ) )
    IARL = ( ( 0 * (IARSmax/(HARS**(1+q))) * RS**(1+q) + pARL * (IARLmax/(HARL**(1+q))) * RL**(1+q) ) / 
    ( 1 + (pARS/(HARS**(1+q))) * RS**(1+q) + (pARL/(HARL**(1+q))) * RL**(1+q) ) )
    # Stage-specific production rates:
    vJ     = betaJ*(IJRS+IJRL) - mJ                 # juvenile production rate
    gammaJ = (vJ - muJ) / (1 - zJA**(1-(muJ/vJ)))   # maturation rate
    vA     = betaA*(IARS+IARL) - mA                 # reproduction rate
    # ODE system:
    dRSdt = deltaRS*(RSmax - RS) - IJRS*J - IARS*A
    dRLdt = deltaRL*(RLmax - RL) - IJRL*J - IARL*A
    dJdt = max(vA,0)*A + vJ*J - max(gammaJ,0)*J - muJ*J
    dAdt = max(gammaJ,0)*J + min(vA,0)*A - muA*A
    return np.array([dRSdt, dRLdt, dJdt, dAdt])

# Initial densities for RS, RL, J, A
X0 = np.array([0.01, 0.01, 0.01, 0.01])

# Time range
t = np.linspace(0,300,1000)

# Solve ODE
X = odeint(model,X0,t)

# Plot results
RS,RL,J,A = np.transpose(X)

plt.figure()
plt.plot(t, RS, 'g-', label='RS', linewidth=1.0)
plt.plot(t, RL, 'g-', label='RL', linewidth=2.5)
plt.legend(loc='upper right')
plt.xlabel('Time (day)')
plt.ylabel('Density (mg/L)')
plt.show()

plt.figure()
plt.plot(t, J, 'k-', label='J', linewidth=1.0)
plt.plot(t, A, 'k-', label='A', linewidth=2.5)
plt.legend(loc='upper right')
plt.xlabel('Time (day)')
plt.ylabel('Density (mg/L)')
plt.show()