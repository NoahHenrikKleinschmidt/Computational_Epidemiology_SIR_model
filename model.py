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

# Alright Battleplan:
# - setup n_subgroups (int) << CHECK
# - setup percentages (cuple of len n) << CHECK
# - setup factors (tuple of len n) << CHECK
# - methods to setup n_subgroups -> init << CHECK
# - method to setup factors << CHECK
# - modify ratio_impact to work with percentages  << CHECK


class SIR:
    """
    This class handles the SIR Model
    """
    def __init__(self, subgroups:int):
        """
        The SIR Model for a given population containing n divergent subgroups. 
        Note that subgroups does not include the "normally" susceptibles (they will be added externally)
        """
        self._pop = 100
        self._infected = 1
        
        self._timestart = 0
        self._timeend = 10
        self._timespace = np.linspace(self._timestart, self._timeend, 10000)

        self._subgroups = subgroups+1

        # setup percentages of the different subgroups (initial all equal)
        self._percentages = [1/self._subgroups for i in range(self._subgroups)]    

        self._initials = (100, 1, 0, 0) # initial values for SIR model (100 people, 1 infected, 0 recovered or dead)

        self._inf_rate = 0.4        # infection rate
        self._rec_rate = 0.7        # recovery rate
        self.__death_rate = 0.01    # death rate
        self._rel_rate = 0.005      # relapsation rate

        # setup the rates factors for the different subgroups
        self._inf_factor = [1 for i in range(subgroups)]
        self._rec_factor = [1 for i in range(subgroups)]
        self._death_factor = [1 for i in range(subgroups)]
        self._rel_factor = [1 for i in range(subgroups)]


    def groups(self):
        """
        Returns the number of subgroups
        """
        return self._subgroups

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
        Set new scaling factors for the different subgroups
        (NOTE: factors are tuples of length subgroups-1, as the 
        first group are the reference group of "normally" susceptibles)
        """

        str_attributes = ["_inf_factor", "_rec_factor", "_death_factor", "_rel_factor"]
        factors = [infection_factor, recovery_factor, death_factor, relapsation_factor]
        
        for factor, attribute in zip(factors, str_attributes):
            
                self._update_factor(attribute, factor)

        return self._inf_factor, self._rec_factor, self._death_factor, self._rel_factor
  

    def percentages(self, p:tuple=None):
        """
        Set the percentages for the different subgroups
        (has to have an entry for each divergent subgroup, but NOT the "normally" susceptibles!)
        """
        if p is not None:
            # convert percentages into list
            p = [p] if isinstance(p, (float, int)) else list(p)
            if len(p) == self._subgroups-1 and sum(p) <= 1:
                p.insert(0, 1-sum(p)) # add normally susceptible percentage at the beginning
                self._percentages = p
            else:
                print(f"Percentages tuple has to have {self._subgroups} entries (got {len(p)}), and sum up to max 1 (current sum = {sum(p)})...")
        else:
            return self._percentages

    def _update_factor(self, attr:str, factor:tuple):
        """
        Update a factor attribute based on its name as as string
        """             
        if factor is None: 
            return
        factor = [factor] if isinstance(factor, (float, int)) else list(factor)
        if len(factor) == self._subgroups-1:
            setattr(self, attr, factor)
          

    def _ratio_impact(self, factor:tuple):
        """
        The Φ impact of highly susceptibles
        """
        impact = self._percentages[0] 
        for p, f in zip(self._percentages[1:], factor):
            impact += f * p
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

if __name__ == '__main__':

    model = SIR(subgroups = 2)
    model.percentages((0.43, 0.22))
    model.factors(
        infection_factor = (3, 12),
        recovery_factor = (12, 0.99),
    )
    print(model.factors())

    t, s = model.solve()
    print(s)