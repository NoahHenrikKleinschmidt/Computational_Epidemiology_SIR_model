"""
----------------------------------------------------------------------------
This module handles interactive charts to visualise results of model.solve()
Note, that these cannot be used to visualse the results of simulate() 
----------------------------------------------------------------------------
"""
import plotly.graph_objects as go
import pandas as pd
import numpy as np 


def LineChart(model):
    """
    Generates a line chart with X-timespan and 
    Y-the proportion of individuals for S, I, R, 
    and D as separate lines.
    """
    fig = go.Figure()

    timespace, solutions = model.timespace(), model.solutions()

    for sol, name in zip(solutions, ["Susceptibles", "Infected", "Recovered", "Dead"]):
        fig.add_trace(
            go.Scatter(
                        x=timespace, y=sol,
                        mode="lines",
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
            title = "Proportion of Individuals per Category"
        )
    )
    return fig

def TimepointBarChart(t, model):
    """
    Generates a bar chart with the distribution of 
    Susceptibles, Infected, Recovered, and Dead at a given timepoint t
    """
    fig = go.Figure()

    timepoints = model.timespace()
    solutions = model.solutions()
    population = sum(model.initials())


    # get index of query timepoint (if not found, check for alternative assignment)
    idx = np.where(timepoints == t)
    if len(idx[0]) == 0: 
        if t > 1:
            idx = t
        else:
            idx = int(len(timepoints) * t)
    else:
        idx = idx[0]

    # timepoint used 
    try: 
        tpoint = round(timepoints[idx],2)
    except: 
        tpoint = t

    names = ["Susceptibles", "Infected", "Recovered", "Dead"]

    df = pd.DataFrame(
                    dict(
                        names = names, 
                        y = [round(i[idx][0], 2) for i in solutions],
                    )
                )       
    fig.add_trace(
        go.Bar(     
            
                    x=df["names"], 
                    y=df["y"],
                    # name=names, 
                    hoverinfo = "y"
                )
            )

        

    fig.update_layout(
        yaxis = dict( range = (0, population)), 
        title = f"Population at t = {tpoint}",
    )
    return fig

def RValueChart(model):
    """
    Generates a line chart of R values over timespan
    """
    fig = go.Figure()
    
    timespace = model.timespace()
    Rvals = model.R()
    
    fig.add_trace(
        go.Scatter(
                    x=timespace, y=Rvals,
                    mode="lines",
                    hoverinfo = "y"
                )
            )
    fig.update_layout(
        title = "R Value",
        xaxis = dict(
            title = "Timespan (e.g. weeks or months)"
        ), 
        # yaxis = dict(
        #     title = "Number of People per Category"
        # )
    )
    return fig

def PhaseChart(model):
    """
    Generates a Phase Chart (scatterplot) with 
    X-Susceptibles and Y-Infectuous over timespan (colorscale)
    """
    fig = go.Figure()
    
    timespace = model.timespace()
    susceptibles, infectuous = model.solutions()[0], model.solutions()[1]
    
    fig.add_trace(
        go.Scatter(
                    x=susceptibles, y=infectuous,
                    mode="markers",
                    hoverinfo = "y+x",
                    marker = dict(
                        cmin = min(timespace), 
                        cmax = max(timespace),
                        color=timespace,
                        colorbar = dict( title = "Timespan \n(e.g. weeks or months)" ),
                        colorscale="Viridis",
                    ),
                )
            )
            
    fig.update_layout(
        title = "Phase Chart",
        xaxis = dict(
            title = "Proportion of Susceptibles"
        ), 
        yaxis = dict(
            title = "Proportion of Infectuous"
        )
    )
    return fig