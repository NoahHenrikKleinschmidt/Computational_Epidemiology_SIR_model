import model
import plotly.graph_objs as go
import streamlit as st




st.title("Impact of Population Heterogeneity for Disease Dynamics")

st.markdown("---")

controls_container = st.container()
controls_panel = controls_container.expander("Parameter Controls")
plot_container = st.container()
line_chart, bar_chart = plot_container.columns((2,1))
bottom_container = st.container()

c1, c2, c3 = controls_panel.columns((1, 1, 1))

# setup parameter controls

c1.markdown("Population Settings \n ---")

high_p = c1.slider("Percentage of highly susceptibles", min_value = 0.001, max_value = 0.999, value = 0.02, step = 0.001)
model.norm_p = 1 - high_p

starting_population = c1.slider("Population size", min_value = 10, max_value = 1000, value = 100)
starting_infectuous = c1.slider("Initial Spreaders", min_value = 0, max_value = 100, value = 1)

c2.markdown("Disease Rates \n ---")

infection_rate = c2.slider("Infection Rate", min_value = 0.01, max_value = 0.15, value = 0.08, step = 0.001, help = "The rate at which susceptibles become infected. Conceptually, the speed-of-spread.")
model.infection_rate = infection_rate

recovery_rate = c2.slider("Recovery Rate", min_value = 0.01, max_value = 2.0, value = 0.9, step = 0.01, help = "The rate at which infectous people recover. The percentage of infectous people who overcame the disease in a given timeframe (e.g. per day). Conceptually, a speed-of-recovery." )
model.recovery_rate = recovery_rate

theta = c2.slider("Death Rate", min_value = 0.0, max_value = 0.05, value = 0.004, step = 0.001, format = "%e", help = "The rate at which infectous people die.\nConceptually, the percentage of infectous people whose disease status changed, and for whom this change meant death (rather than recovery).")
model.theta = theta

relapsation_rate = c2.slider("Relapsation Rate", min_value = 0.000, max_value = 0.9, value = 0.01, step = 0.001, format = "%e", help = "The rate at which recovered relaps to susceptibles.\nConceptually similar to the percentage of people who loose their immunity every day.")
model.relapsation_rate = relapsation_rate

c3.markdown("Highly Susceptible's \n --- ")


k = c3.slider("Infection Risk", min_value = 0.01, max_value = 10.0, value = 1.5, step = 0.1, help = "Elevated infection risk. Higher values mean higher chance of infection.")
model.k = k

j = c3.slider("Recovery Delay", min_value = 0.01, max_value = 5.0, value = 1.1, step = 0.01, help = "Delay in recovery. Higher values mean longer time to recovery.")
model.j = 1/j

q = c3.slider("Death Risk", min_value = 0.01, max_value = 10.0, value = 1.1, step = 0.1, help = "Elevated death risk. Higher values mean increased risk of dying (instead of recovering).")
model.q = q

h = c3.slider("Relapsation Risk", min_value = 0.01, max_value = 10.0, value = 1.3, step = 0.1, help = "Elevated relapsation rate. Higher values mean faster loss of immunity.")
model.h = h


initials = ( \
    starting_population,     # S
    starting_infectuous,      # I
    0,      # R
    0,      # D 
)

population = starting_population + starting_infectuous

endpoint = bottom_container.slider("End Timepoint", min_value = 2, max_value = 100, step = 1, value = 10)

timespace, solutions = model.solve(initials, end = endpoint)
linechart = model.LineChart(timespace, solutions)
line_chart.plotly_chart(linechart, use_container_width=True)
barchart = model.TimepointBarChart(endpoint, timespace, solutions, population)
bar_chart.plotly_chart(barchart, use_container_width=True)