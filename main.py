import model
import plotly.graph_objs as go
import streamlit as st


# initialise the model
sir = model.SIR()

# setup app 

st.set_page_config(
     page_title="Disease Dynamics SIR Modeller",
     page_icon="📈",
     layout="wide",
    #  menu_items={
    #      'To the Repo': 'https://github.com/NoahHenrikKleinschmidt/Computational_Epidemiology_SIR_model',
         
    #  }
)

st.sidebar.markdown("""
The underlying system of equations can be found at the repo of this app on [Github](https://github.com/NoahHenrikKleinschmidt/Computational_Epidemiology_SIR_model)
""")

st.sidebar.markdown("## Set slider limits \n ---")

# set slider limits
max_population = st.sidebar.number_input("Population Size", value = 1000)

st.title("Impact of Population Heterogeneity for Disease Dynamics")

st.markdown("---")

# setup layout
controls_container = st.container()
controls_panel = controls_container.expander("Parameter Controls")
plot_container = st.container()
line_chart, bar_chart = plot_container.columns((2,1))
bottom_container = st.container()
c1, c2, c3 = controls_panel.columns((1, 1, 1))


# setup parameter controls
c1.markdown("Population Settings \n ---")

high_p = c1.slider("Percentage of highly susceptibles", min_value = 0.001, max_value = 0.999, value = 0.02, step = 0.001)
starting_population = c1.slider("Population size", min_value = 10, max_value = max_population, value = 100)
starting_infectuous = c1.slider("Initial Spreaders", min_value = 0, max_value = 100, value = 1)

# add values to the sir model
sir.percentage(high_p)
sir.initials( values = (starting_population, starting_infectuous, 0, 0) )


c2.markdown("Disease Rates \n ---")

infection_rate = c2.slider("Infection Rate", min_value = 0.01, max_value = 0.15, value = 0.08, step = 0.001, help = "The rate at which susceptibles become infected. Conceptually, the speed-of-spread.")
recovery_rate = c2.slider("Recovery Rate", min_value = 0.01, max_value = 2.0, value = 0.9, step = 0.01, help = "The rate at which infectous people recover. The percentage of infectous people who overcame the disease in a given timeframe (e.g. per day). Conceptually, a speed-of-recovery." )
theta = c2.slider("Death Rate", min_value = 0.0, max_value = 0.05, value = 0.004, step = 0.001, format = "%e", help = "The rate at which infectous people die.\nConceptually, the percentage of infectous people whose disease status changed, and for whom this change meant death (rather than recovery).")
relapsation_rate = c2.slider("Relapsation Rate", min_value = 0.000, max_value = 0.9, value = 0.01, step = 0.001, format = "%e", help = "The rate at which recovered relaps to susceptibles.\nConceptually similar to the percentage of people who loose their immunity every day.")

# add values to the sir model
sir.rates(
    infection_rate = infection_rate, 
    recovery_rate = recovery_rate,
    death_rate = theta,
    relapsation_rate = relapsation_rate
)


c3.markdown("Highly Susceptible's \n --- ")

infection_factor = c3.slider("Infection Risk", min_value = 0.01, max_value = 10.0, value = 1.5, step = 0.1, help = "Elevated infection risk. Higher values mean higher chance of infection.")
recovery_factor = c3.slider("Recovery Delay", min_value = 0.01, max_value = 5.0, value = 1.1, step = 0.01, help = "Delay in recovery. Higher values mean longer time to recovery.")
death_factor = c3.slider("Death Risk", min_value = 0.01, max_value = 10.0, value = 1.1, step = 0.1, help = "Elevated death risk. Higher values mean increased risk of dying (instead of recovering).")
relapsation_factor = c3.slider("Relapsation Risk", min_value = 0.01, max_value = 10.0, value = 1.3, step = 0.1, help = "Elevated relapsation rate. Higher values mean faster loss of immunity.")

# add values to the sir model
sir.factors(
    infection_factor = infection_factor, 
    recovery_factor = 1/recovery_factor,
    death_factor = death_factor,
    relapsation_factor = relapsation_factor
)


# setup general properties of the simulation

population = starting_population + starting_infectuous
endpoint = bottom_container.slider("End Timepoint", min_value = 2, max_value = 100, step = 1, value = 10)

# add values to the sir model
sir.set_timespace(stop = endpoint)



# solve the SIR Model
timespace, solutions = sir.solve()


# visualise Results

linechart = model.LineChart(timespace, solutions)
line_chart.plotly_chart(linechart, use_container_width=True)
barchart = model.TimepointBarChart(endpoint, timespace, solutions, population)
bar_chart.plotly_chart(barchart, use_container_width=True)