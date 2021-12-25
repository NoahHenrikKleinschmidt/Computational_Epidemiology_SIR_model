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

---------------------------------------------------------------------------
Usage:
---------------------------------------------------------------------------
# setup 
model = SIR(subgroups = 2) # 2 divergent subgroups

# defining initial parameters
model.rates( infection_rate = 0.04 )
 ...

# solving the model equations for a fix percentage ratio of subgroups
timespace, solutions = model.solve()

# visualising results
fig = charts.LineChart(timespace, solutions)
fig.show()


# or simulate changin percentage ratios of subgroups
# (Note: this corresponds to varying initial percentages of subgroups, 
#        NOT of dynamically incrementing or decreasing their 
#        precentages...)

# (start_percentage, end_percentage) for each subgroup
percentage_array = [(0.1, 0.4), (0.3, 0.5)] 

# generates a matplotlib figure (also models development of R value), 
# divides start-end into interval of 10-steps
simulate(model, p = percentage_array, steps = 10, model_R = True) 

---------------------------------------------------------------------------
"""

from scipy.integrate import solve_ivp # numerical ODE solver from scipy
import numpy as np 
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib.lines import Line2D
import pandas as pd
import interactive_charts as charts

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
        
        self._timestart = 0
        self._timeend = 20
        self._timespace = np.linspace(self._timestart, self._timeend, 10000)

        self._solutions = None

        self._subgroups = subgroups+1

        # setup percentages of the different subgroups (initial all equal)
        self._percentages = [1/self._subgroups for i in range(self._subgroups)]    

        self._initials = (0.98, 0.02, 0, 0) # initial values for SIR model (98% people susceptible, 2% infected, 0 recovered or dead)

        self._inf_rate = 0.4        # infection rate
        self._rec_rate = 0.6        # recovery rate
        self.__death_rate = 0.05    # death rate
        self._rel_rate = 0.008      # relapsation rate

        # setup the rates factors for the different subgroups
        self._inf_factor = [1 for i in range(subgroups)]
        self._rec_factor = [1 for i in range(subgroups)]
        self._death_factor = [1 for i in range(subgroups)]
        self._rel_factor = [1 for i in range(subgroups)]

    def solutions(self):
        """
        Returns the solutions of self.solve()
        """
        return self._solutions

    def groups(self):
        """
        Returns the number of subgroups
        """
        return self._subgroups

    def model(self, t, initials):
        """
        The SIRD Model
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
        self._solutions = solutions
        return self._timespace, solutions

    def R0(self):
        """
        Returns R0 as weighted (infection_rate + relapsation_rate) / (recovery_rate + death_rate)
        """
        upper = self._infection_rate() + self._relapsation_rate()
        lower = self._recovery_rate() + self._death_rate()
        R0 = upper / lower
        return R0

    def R(self, S:(float or np.array) = None):
        """
        Returns the R value either for a single timepoint 
        (if S = float) or an entire timespan (if S = np.array), 
        if S is None, precomputed solutions will be used...
        """
        if S is None: 
            S = self.solutions()[0]
        elif isinstance(S, (float)):
            if S > 1: 
                S = S / sum(self._initials) # transform into proportion within the population
        else: 
            if any([i > 1 for i in list(S)]):
                S = S / sum(self._initials)

        R = S * self.R0()
        return R

    def initials(self, values=None):
        """
        Setup initial values for the model
        Default: (0.98, 0.02, 0, 0)
        """
        if values is None: 
            return self._initials
        else: 
            # convert to proprtion
            if any([i > 1 for i in values]):
                values = tuple([i / sum(values) for i in values])
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

def simulate(model, p:np.array, steps = 10, model_R = False, summary = False, ax = None, show=False, **kwargs):
    """
    Simulate the effect of changing population percentages 
    within a model, where p is an array of starting and end percentages of each subgroup as tuple.
    """

    # get some kwargs and setup default values
    show_legend = kwargs.pop("show_legend", True)
    colors = kwargs.pop("colors", ["mediumblue", "darkred", "darkcyan", "black"])
    yscale = kwargs.pop("yscale", "linear")
    show_threshold = kwargs.pop("threshold", True)
    first_alpha = kwargs.pop("fadeout", 0.4)
    show_plot_titles = kwargs.pop("plot_titles", True)

    # prepare plt subplots for plotting
    ax, Rax, sum_ax = _prepare_axes(model_R, summary, ax, kwargs)

    # summary endpoint statistics dict    
    end_stats = {"Susceptibles" : [], "Infectuous" : [], "Recovered" : [], "Deceased" : []}

    # generate percentages for steps and corresponding fadeout alphas
    percents = [np.linspace(start = i[0], stop = i[1], num = steps) for i in list(p)]
    percents = zip(*percents) # transpose to get the ith percentage for each subgroup together...
    alphas = [1] + list(np.linspace(first_alpha, 0.01, num = steps))


    # make the main linechart (+ R chart if model_R)
    # yeah, the arguments are really long...
    _dynamics_linechart(model, model_R, summary, ax, kwargs, show_legend, colors, yscale, show_threshold, show_plot_titles, Rax, end_stats, percents, alphas)

    # make additional summary chart
    if summary:
        _summary_barchart(sum_ax, end_stats, yscale, show_plot_titles)

    if show: 
        plt.tight_layout()
        plt.show()



def _dynamics_linechart(model, model_R, summary, ax, kwargs, show_legend, colors, yscale, show_threshold, show_plot_titles, Rax, end_stats, percents, alphas):
    """
    Draw the main disease dynamics line chart
    """
    for percent, alpha in zip(percents, alphas):
        model.percentages(percent)

        t, sol = model.solve()

        for s, c in zip(sol, colors): 
            ax.plot(t, s, c = c, alpha = alpha, **kwargs)

        if summary: 
            end_stats["Susceptibles"].append(sol[0][-1])
            end_stats["Infectuous"].append(sol[1][-1])
            end_stats["Recovered"].append(sol[2][-1])
            end_stats["Deceased"].append(sol[3][-1])

        if model_R:
            Rax.plot(t, model.R(sol[0]), alpha = alpha, c = "darkslategray")


    # graph formatting
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_yscale(yscale) # log scale makes seeing the different lines more easy
    ax.set(
        xlabel = "Timespan (e.g. weeks or months)", 
        title = "Disease Dynamics" if show_plot_titles else None, 
        ylabel = "Proportion of \n individuals per Category",
    )
    
    if show_threshold:
        _immunity_thershold_line(model, ax)

    if show_legend:
        ax.legend(
            handles = [
                        Line2D([0], [0], color=colors[0], lw=2, label="Susceptibles"),
                        Line2D([0], [0], color=colors[1], lw=2, label="Infectuous"),
                        Line2D([0], [0], color=colors[2], lw=2, label="Recovered"),
                        Line2D([0], [0], color=colors[3], lw=2, label="Deceased"),
                    ],
            bbox_to_anchor = (1, 1), frameon = False
        )

    if model_R: # some more reformatting of the model_R chart
        _Rvalue_chart_reformat(Rax, model, yscale, show_threshold, show_plot_titles)

def _immunity_thershold_line(model, ax):
    """
    Draw a hashed herd immunity line
    """
    # population threshold for herd immunity (HIT) (for susceptibles)
    ax.hlines(
            y = 1 / model.R0() , 
            xmin = min(model.timespace()), 
            xmax = max(model.timespace()), 
            color = "dimgray", linestyle = "--", linewidth = 2,
        )


def _Rvalue_chart_reformat(ax, model, yscale, show_threshold, show_plot_titles):
    """
    Generate a linechart showing the dynamics of the R value
    """
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_yscale(yscale) # log scale makes seeing the different lines more easy
    ax.set(
            xlabel = "Timespan (e.g. weeks or months)", 
            title = "R Value Dynamics" if show_plot_titles else None, 
            ylabel = "R value"
        )

    if show_threshold:
        # R == 1
        ax.hlines(
                    y = 1,  
                    xmin = min(model.timespace()), 
                    xmax = max(model.timespace()), 
                    color = "dimgray", linestyle = "--", linewidth = 2,
                )

def _summary_barchart(ax, end_stats, yscale, show_plot_titles):
    """
    Generate an endpoint population summary barchart...
    """
    df = pd.DataFrame(end_stats)
    df = df.transpose()
    df.plot.bar(    
                        ax = ax, 
                        colormap = "Blues_r", 
                        legend = False, rot = 0,
                        edgecolor = "black",
                        alpha = 0.8, 
                        )
    # axis formatting
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_yscale(yscale) # log scale makes seeing the different lines more easy
    ax.set(
            xlabel = "", 
            title = "Population at Endpoint of Simulation" if show_plot_titles else None,
            ylabel = "Proportion of \n individuals per Category",
        )


def _prepare_axes(model_R, summary, ax, kwargs):
    """
    Extracts the subplots for plotting or generates a new 
    figure with the appropriate number of subplots
    """
    # check what additional plots must be made
    # and see if axs have already been passed for those
    make_extra_plots = False
    try:
        if model_R:
            Rax = ax[1]
        if summary:
            sum_ax = ax[2] if model_R else ax[1]
        if summary or model_R:
            ax = ax[0]

    except:
        make_extra_plots = True

    # if additional axes are required, make a new figure...
    # this will ignore any submitted axes!

    if ax is None or make_extra_plots: 
        # get number of required plots
        subplots = 1 + summary + model_R

        figsize = kwargs.pop("figsize", None)
        fig, axs = plt.subplots(subplots, figsize = figsize)

        ax = axs[0] if subplots > 1 else axs
        Rax = axs[1] if model_R else None
        if summary: 
            sum_ax = axs[2] if model_R else axs[1]
        else:
            sum_ax = None

    return ax,Rax,sum_ax


if __name__ == "__main__":
    model = SIR(subgroups = 2)
    model.initials(
        (98, 12, 0, 0)
    )
    model.set_timespace(stop = 100)

    model.rates( infection_rate = 0.5, recovery_rate = 0.5 )

    model.factors( 
        infection_factor    =       (20, 0.4),
        recovery_factor     =       (0.5, 2), 
        death_factor        =       (4, 0.1),
        relapsation_factor  =       (2.5, 1),
    )

    # fig, axs = plt.subplots(ncols = 2)

    # ax1 = simulate(
    #                 model, 
    #                 p = np.array([(0.10, 0.75), (0.75, 0.25)]), 
    #                 ax = axs[0], R_ax = axs[1],
    #                 model_R = True,
    #                 # yscale = "log", 
    #             )

    # plt.tight_layout()
    # plt.show()

    model.solve()
    t = .15
    rval_chart = charts.LineChart(model)
    rval_chart.show()
