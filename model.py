"""
This module specifies the SIR Model for this project
Parameters of this model are: 

---------------------------------------------------------------------------
norm_p      The percentage of normally susceptibles within the population
*_rate      A transition rate for normally susceptibles,
            highly susceptibles will be factor x *_rate 
k           k-fold increased infection rate for highly susceptibles (vs. 
            normally susceptibles)
theta       percentage of I→D from I→(R'+D)
j           k-fold decreased recovery rate for highly susceptibles
h           h-fold changed relapsation rate for highly susceptibles
q           q-fold increased death rate for highly susceptibles
---------------------------------------------------------------------------
"""

from scipy.integrate import solve_ivp # numerical ODE solver from scipy
import numpy as np 
import matplotlib.pyplot as plt
import plotly.graph_objects as go




def _infection_rate(infection_rate, norm_p, k):
    """
    Returns the cumulated rate at which susceptibles transition to infected
    """
    rate = 1 - norm_p
    rate = norm_p + k * rate
    rate = infection_rate * rate
    return rate

def _relapsation_rate(relapsation_rate, h):
    """
    Returns the rate at which Recovered people relapse
    """
    rate = relapsation_rate * ( 1 + h )
    return rate

def _recovery_rate(recovery_rate, theta, j):
    """
    Returns the rate at which infected transition to recovered
    """
    rate = (1 - theta) * recovery_rate
    rate = ( 1 + j ) * rate
    return rate

def _death_rate(theta, q):
    """
    Returns the rate at which infected transition to dead
    """
    rate = (1 + q)*theta
    return rate

# setup the surrounding parameters of the SIRD Model

norm_p = 0.95
infection_rate = 0.7    # infection rate
k = 1.6    # k = elevated infection rate for highly susceptibles
recovery_rate = 0.5    # recovery rate
j = 0.8    # j = reduced recovery rate for highly susceptibles
theta = 0.2    # death rate (percentage of recovery rate that are actually dying not recovering as such)
q = 1.2    # q = elevated death rate for highly susceptibles
relapsation_rate = 0.1    # relapsation rate
h = 1.3    # h = elevated relapsation rate for highly susceptibles

# setup the system of equations
def SIRD(
        t, 
        initials:tuple
    ):
    
    # setup initial conditions
    S = initials[0]
    I = initials[1]
    R = initials[2]
    D = initials[3]



    # setup the diff equations
    dS_dt = - _infection_rate(infection_rate, norm_p, k) * S * I \
            + _relapsation_rate(relapsation_rate, h) * R


    dI_dt = _infection_rate(infection_rate, norm_p, k) * S * I \
            - _recovery_rate(recovery_rate, theta, j) * I \
            - _death_rate(theta, q) * I
    
    dR_dt = _recovery_rate(recovery_rate, theta, j) * I \
            - _relapsation_rate(relapsation_rate, h) * R
    
    dD_dt = _death_rate(theta, q) * I

    return [dS_dt, dI_dt, dR_dt, dD_dt]

# now setup simulation
startpoint = 0 
endpoint = 100
steps = 10000
timespace = np.linspace(startpoint, endpoint, steps)

initials = ( \
    99995,  # S
    5,      # I
    0,      # R
    0,      # D 
)

# solve SIRD_system
SIRD_system = solve_ivp(
                    SIRD, 
                    (startpoint, endpoint), 
                    initials, 
                    method="RK45", 
                    dense_output=True
                )

# get and transpose results
solution = SIRD_system.sol(timespace)
solutions = [i.T for i in solution]

# visualise results

def LineChart():
    fig = go.Figure()

    for sol, name in zip(solutions, ["S", "I", "R", "D"]):
        fig.add_trace(
            go.Scatter(
                        x=timespace, y=sol,
                        mode='lines',
                        name=name
                    )
                )

    return fig

fig = LineChart()
fig.show()