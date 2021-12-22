"""
This module specifies the SIR Model for this project

---------------------------------------------------------------------------
Parameters of this model are: 
---------------------------------------------------------------------------
percentage  The percentage of highly susceptibles within the population 
            (stored as 1-p, since the model tracks the percentage of 
            normally susceptibles...)
*_rate      A transition rate for normally susceptibles,
            highly susceptibles will be factor x *_rate 
*_factor    The scalar factor for highly susceptibles
initials    Initial values for population
            Default: 100 susceptibles, 1 infectuous, 0 dead, 0 recovered
---------------------------------------------------------------------------

-----------------------------------------
Usage:
-----------------------------------------
# setup 
model = SIR()

# defining initial parameters
model.rates( infection_rate = 0.04 )
 ...

# solving the model equations
timespace, solutions = model.solve()

# visualising results
fig = LineChart(timespace, solutions)
fig.show()
-----------------------------------------
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