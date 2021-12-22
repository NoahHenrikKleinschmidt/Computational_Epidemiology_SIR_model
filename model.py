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
import pandas as pd


def ratio_impact(norm_p, factor):
    """
    Returns a_n + factor*a_h, 
    the Φ impact of highly susceptibles
    """
    impact = norm_p + factor * (1 - norm_p)
    return impact

def _infection_rate(infection_rate, norm_p, k):
    """
    Returns the cumulated rate at which susceptibles transition to infected
    """
    rate = ratio_impact(norm_p, k)
    rate = rate * infection_rate
    return rate

def _relapsation_rate(relapsation_rate, norm_p, h):
    """
    Returns the rate at which Recovered people relapse
    """
    rate = ratio_impact(norm_p, h)
    rate = relapsation_rate * rate
    return rate

def _recovery_rate(recovery_rate, norm_p, j):
    """
    Returns the rate at which infected transition to recovered
    """
    rate = ratio_impact(norm_p, j)
    rate = rate * recovery_rate
    return rate

def _death_rate(theta, norm_p, q):
    """
    Returns the rate at which infected transition to dead
    """
    rate = ratio_impact(norm_p, q)
    rate = theta * rate
    return rate

# setup the surrounding parameters of the SIRD Model

norm_p = 0.95
infection_rate = 0.4    # infection rate
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
            + _relapsation_rate(relapsation_rate, norm_p, h) * R


    dI_dt = _infection_rate(infection_rate, norm_p, k) * S * I \
            - _recovery_rate(recovery_rate, norm_p, j) * I \
            - _death_rate(theta, norm_p, q) * I
    
    dR_dt = _recovery_rate(recovery_rate, norm_p, j) * I \
            - _relapsation_rate(relapsation_rate, norm_p, h) * R
    
    dD_dt = _death_rate(theta, norm_p, q) * I

    return [dS_dt, dI_dt, dR_dt, dD_dt]


# initials = ( \
#     99995,  # S
#     5,      # I
#     0,      # R
#     0,      # D 
# )

def solve(initials, start = 0, end = 100, steps = 10000):
    """
    Solve the equations based on initial parameters...
    """
    # now setup simulation
    timespace = np.linspace(start, end, steps)

    # solve SIRD_system
    SIRD_system = solve_ivp(
                        SIRD, 
                        (start, end), 
                        initials, 
                        method="RK45", 
                        dense_output=True
                    )

    # get and transpose results
    solution = SIRD_system.sol(timespace)
    solutions = [i.T for i in solution]
    return timespace, solutions

# visualise results

def LineChart(timespace, solutions):
    fig = go.Figure()

    for sol, name in zip(solutions, ["Susceptibles", "Infected", "Recovered", "Dead"]):
        fig.add_trace(
            go.Scatter(
                        x=timespace, y=sol,
                        mode='lines',
                        name=name, 
                        hoverinfo = "y+name"
                    )
                )
    fig.update_layout(
        title = "Disease Dynamics",
        xaxis = dict(
            title = "Timespan (e.g. weeks or months)"
        ), 
        yaxis = dict(
            title = "Number of People per Category"
        )
    )
    return fig

def TimepointBarChart(t, timepoints, solutions, population):
    fig = go.Figure()

    idx = np.where(timepoints == t)

    names = ["Susceptibles", "Infected", "Recovered", "Dead"]

    df = pd.DataFrame(
                    dict(
                        names = names, 
                        y = [round(i[idx][0]) for i in solutions],
                    )
                )
    fig.add_trace(
        go.Bar(     
            
                    x=df["names"], 
                    y=df["y"],
                    # name=names, 
                    # hoverinfo = "y+name"
                )
            )
    fig.update_layout(
        yaxis = dict( range = (0, population)), 
        title = f"Population at t = {t}",
    )
    return fig

# fig = LineChart()
# fig.show()