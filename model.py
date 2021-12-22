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


class SIR:
    """
    This class handles the SIR Model
    """
    def __init__(self):
        self._pop = 100
        self._infected = 1
        
        self._timestart = 0
        self._timeend = 10
        self._timespace = np.linspace(self._timestart, self._timeend, 10000)

        self._percent_norms = 0.95      # percentage of "normally" susceptibles
        self._initials = (100, 1, 0, 0) # initial values for SIR model (100 people, 1 infected, 0 recovered or dead)

        self._inf_rate = 0.4        # infection rate
        self._rec_rate = 0.7        # recovery rate
        self.__death_rate = 0.01    # death rate
        self._rel_rate = 0.005      # relapsation rate

        self._inf_factor = 1.5
        self._rec_factor = 0.6
        self._death_factor = 1.2
        self._rel_factor = 1.4

    def model(self, t, initials):
        """
        The SIR Model
        """
        S = initials[0]
        I = initials[1]
        R = initials[2]
        D = initials[3]

        # setup the diff equations
        dS_dt = - self._infection_rate() * S * I \
                + self._relapsation_rate() * R


        dI_dt = self._infection_rate() * S * I \
                - self._recovery_rate() * I \
                - self._death_rate() * I
    
        dR_dt = self._recovery_rate() * I \
                - self._relapsation_rate() * R
    
        dD_dt = self._death_rate() * I

        return [dS_dt, dI_dt, dR_dt, dD_dt]


    def solve(self):
        """
        Solve the equations based on initial parameters...
        """
        
        # solve SIRD_system
        SIRD_system = solve_ivp(
                            self.model, 
                            (self._timestart, self._timeend), 
                            self._initials, 
                            method="RK45", 
                            dense_output=True
                        )

        # get and transpose results
        solution = SIRD_system.sol(self._timespace)
        solutions = [i.T for i in solution]
        return self._timespace, solutions

    def initials(self, values=None):
        """
        Setup initial values for the model
        Default: (100, 1, 0, 0)
        """
        if values is None: 
            return self._initials
        else: 
            self._initials = values

    def set_timespace(self, start=0, stop=10, step=10000):
        """
        Set a new timespace
        """
        self._timestart, self._timeend = start, stop
        self._timespace = np.linspace(self._timestart, self._timeend, step)

    def timespace(self):
        """
        Returns the timespace
        """
        return self._timespace

    def rates(self, infection_rate=None, recovery_rate=None, death_rate=None, relapsation_rate=None):
        """
        Set new rates (or get current ones)
        """
        if infection_rate is not None: self._inf_rate = infection_rate
        if recovery_rate is not None: self._rec_rate = recovery_rate
        if death_rate is not None: self.__death_rate = death_rate
        if relapsation_rate is not None: self._rel_rate = relapsation_rate
        
        return self._inf_rate, self._rec_rate, self.__death_rate, self._rel_rate

    def factors(self, infection_factor=None, recovery_factor=None, death_factor=None, relapsation_factor=None):
        """
        Set new scaling factors for "highly" susceptibles
        """
        if infection_factor is not None: self._inf_factor = infection_factor
        if recovery_factor is not None: self._rec_factor = recovery_factor
        if death_factor is not None: self._death_factor = death_factor
        if relapsation_factor is not None: self._rel_factor = relapsation_factor

        return self._inf_factor, self._rec_factor, self._death_factor, self._rel_factor

    def percentage(self, p=None):
        """
        Set the percentage of "highly" susceptibles in the population
        """
        if p is not None:
            if p > 1: 
                print(f"p must be <= 1! Received: {p} ...")
                return
            self._percent_norms = 1 - p
        else:
            return 1 - self._percent_norms

    def _ratio_impact(self, factor):
        """
        Returns a_n + factor*a_h, 
        -> i.e. the Φ impact of highly susceptibles
        """
        impact = self._percent_norms + factor * (1 - self._percent_norms)
        return impact

    def _infection_rate(self):
        """
        Returns the cumulated rate at which susceptibles transition to infected
        """
        rate = self._ratio_impact(self._inf_factor)
        rate = self._inf_rate * rate
        return rate

    def _relapsation_rate(self):
        """
        Returns the rate at which Recovered people relapse
        """
        rate = self._ratio_impact(self._rel_factor)
        rate = self._rel_rate * rate
        return rate

    def _recovery_rate(self):
        """
        Returns the rate at which infected transition to recovered
        """
        rate = self._ratio_impact(self._rec_factor)
        rate = self._rec_rate * rate
        return rate

    def _death_rate(self):
        """
        Returns the rate at which infected transition to dead
        """
        rate = self._ratio_impact(self._death_factor)
        rate = self.__death_rate * rate
        return rate

# # setup the surrounding parameters of the SIRD Model

# norm_p = 0.95
# infection_rate = 0.4    # infection rate
# k = 1.6    # k = elevated infection rate for highly susceptibles
# recovery_rate = 0.5    # recovery rate
# j = 0.8    # j = reduced recovery rate for highly susceptibles
# theta = 0.2    # death rate (percentage of recovery rate that are actually dying not recovering as such)
# q = 1.2    # q = elevated death rate for highly susceptibles
# relapsation_rate = 0.1    # relapsation rate
# h = 1.3    # h = elevated relapsation rate for highly susceptibles
        
# # setup the system of equations
# def SIRD(
#         t, 
#         initials:tuple
#     ):
    
#     # setup initial conditions
#     S = initials[0]
#     I = initials[1]
#     R = initials[2]
#     D = initials[3]
    
#     # setup the diff equations
#     dS_dt = - _infection_rate(infection_rate, norm_p, k) * S * I \
#             + _relapsation_rate(relapsation_rate, norm_p, h) * R


#     dI_dt = _infection_rate(infection_rate, norm_p, k) * S * I \
#             - _recovery_rate(recovery_rate, norm_p, j) * I \
#             - _death_rate(theta, norm_p, q) * I
    
#     dR_dt = _recovery_rate(recovery_rate, norm_p, j) * I \
#             - _relapsation_rate(relapsation_rate, norm_p, h) * R
    
#     dD_dt = _death_rate(theta, norm_p, q) * I

#     return [dS_dt, dI_dt, dR_dt, dD_dt]


# # initials = ( \
# #     99995,  # S
# #     5,      # I
# #     0,      # R
# #     0,      # D 
# # )

# def solve(initials, start = 0, end = 100, steps = 10000):
#     """
#     Solve the equations based on initial parameters...
#     """
#     # now setup simulation
#     timespace = np.linspace(start, end, steps)

#     # solve SIRD_system
#     SIRD_system = solve_ivp(
#                         SIRD, 
#                         (start, end), 
#                         initials, 
#                         method="RK45", 
#                         dense_output=True
#                     )

#     # get and transpose results
#     solution = SIRD_system.sol(timespace)
#     solutions = [i.T for i in solution]
#     return timespace, solutions

# # visualise results

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