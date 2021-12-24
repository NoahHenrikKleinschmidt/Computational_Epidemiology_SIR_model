import model
import interactive_charts as charts
import plotly.graph_objs as go
import streamlit as st




# setup app 

st.set_page_config(
     page_title="Disease Dynamics SIR Modeller",
     page_icon="📈",
     layout="wide",
     initial_sidebar_state = "collapsed",
    #  menu_items={
    #      'To the Repo': 'https://github.com/NoahHenrikKleinschmidt/Computational_Epidemiology_SIR_model',
         
    #  }
)

st.sidebar.markdown("""
The underlying system of equations can be found at the project repo of this app on [Gitlab](https://gitlab.com/NoahHenrikKleinschmidt/Computation_Epidemiology_HS21).
The repo of this app itself is hosted on [Github](https://github.com/NoahHenrikKleinschmidt/Computational_Epidemiology_SIR_model).
""")

st.sidebar.markdown("## Set slider limits \n ---")

st.sidebar.markdown("#### Transition Rates")
# set slider limits
max_population = st.sidebar.number_input("Maximum Population Size", value = 1000)
max_spreaders = st.sidebar.number_input("Maximum Initial Spreaders", value = 100)
max_infection_rate = st.sidebar.number_input("Maximum Infection Rate", value = 5.0)
max_recovery_rate = st.sidebar.number_input("Maximum Recovery Rate", value = 3.0)
max_death_rate = st.sidebar.number_input("Maximum Death Rate", value = 1.0)
max_relapsation_rate = st.sidebar.number_input("Maximum Relapsation Rate", value = 1.0)

st.sidebar.markdown("#### \"Highly\" Susceptible's scalar factors")
max_infection_factor = st.sidebar.number_input("Maximum Infection Factor", value = 10.0)
max_recovery_factor = st.sidebar.number_input("Maximum Recovery Factor", value = 5.0)
max_death_factor = st.sidebar.number_input("Maximum Death Factor", value = 10.0)
max_relapsation_factor = st.sidebar.number_input("Maximum Relapsation Factor", value = 10.0)


# Battleplan: 
# - Add a slider for n subgroups
# - Make a function to add n control columns for factors
# - store that information in the session_state of streamlit...
# - store assigned factors for each subgroup control column in a tuple that can be assigned to the SIR model based on index of the column...


st.title("Impact of Population Heterogeneity for Disease Dynamics")

st.markdown("---")

# setup layout
controls_container = st.container()
setup_panel = controls_container.expander("Basic SIRD Controls")
controls_panel = controls_container.expander("Subgroup Controls")
plot_container = st.container()
line_controls, line_chart, bar_chart = plot_container.columns((0.2, 2,1))
bottom_container = st.container()

c1, c2  = setup_panel.columns((0.75, 1))

# needs to be redefined for n columns...
# c1, c2, c3 = controls_panel.columns((1, 1, 1))


# setup parameter controls
# controls_container.markdown("Basic SIR Settings \n ---")
n_subgroups = c1.number_input("Number of divergent subgroups", min_value = 1, value = 1, help = "Number of divergent subgroups that shall be created in the population")
st.session_state.n_subgroups = n_subgroups # this is the session state thing that makes the app remember attributes ...

starting_population = c1.slider("Population size", min_value = 2, max_value = int(max_population), value = 100)
starting_infectuous = c1.slider("Initial Spreaders", min_value = 1, max_value = int(max_spreaders), value = 1)


c2.markdown("Disease Rates \n ---")

infection_rate = c2.slider("Infection Rate", min_value = 0.01, max_value = max_infection_rate, value = 1.4, step = 0.01, help = "The rate at which susceptibles become infected. Conceptually, the speed-of-spread.")
recovery_rate = c2.slider("Recovery Rate", min_value = 0.01, max_value = max_recovery_rate, value = 0.9, step = 0.01, help = "The rate at which infectous people recover. The percentage of infectous people who overcame the disease in a given timeframe (e.g. per day). Conceptually, a speed-of-recovery." )
theta = c2.slider("Death Rate", min_value = 0.0, max_value = max_death_rate, value = 0.005, step = 0.001, format = "%e", help = "The rate at which infectous people die.\nConceptually, the percentage of infectous people whose disease status changed, and for whom this change meant death (rather than recovery).")
relapsation_rate = c2.slider("Relapsation Rate", min_value = 0.000, max_value = max_relapsation_rate, value = 0.01, step = 0.001, format = "%e", help = "The rate at which recovered relaps to susceptibles.\nConceptually similar to the percentage of people who loose their immunity every day.")



# generate controls for percentages  
# @st.cache(suppress_st_warning=True) # pretty cool for performance improvement... ^^
def _setup_subgroup_controls(n):
    """
    Make factor control columns for each subgroup
    """

    # setup factors array for all subgroups
    st.session_state.factors_array = [[1 for i in range(4)] for j in range(n)]

    # setup percentages array for all subgroups
    st.session_state.percentages = [1/st.session_state.n_subgroups for i in range(n)]

    # setup column containers for the sliders of each subgroup
    subgroup_control_columns = controls_panel.columns(tuple([1 for i in range(n)]))
    
    # now setup the actual columns
    for c, group in zip(subgroup_control_columns, range(st.session_state.n_subgroups)):
        # make a new control column
        c.markdown(f"#### Subgroup {group+1}\n ---")
        st.session_state.percentages[group] = c.slider("Percentage of the subgroup within the population", key = f"g{group}_Percentage of the subgroup within the population", min_value = 0., max_value = 1., step  = 0.1)
        st.session_state.factors_array[group][0] = c.slider("Infection Risk", key = f"g{group}_Infection Risk", min_value = 0.01, max_value = max_infection_factor, value = 1.5, step = 0.1, help = "Elevated infection risk. Higher values mean higher chance of infection.")
        st.session_state.factors_array[group][1] = c.slider("Recovery Delay", key = f"g{group}_Recovery Delay", min_value = 0.01, max_value = max_recovery_factor, value = 1.1, step = 0.01, help = "Delay in recovery. Higher values mean longer time to recovery.")
        st.session_state.factors_array[group][2] = c.slider("Death Risk", key = f"g{group}_Death Risk", min_value = 0.01, max_value = max_death_factor, value = 1.1, step = 0.1, help = "Elevated death risk. Higher values mean increased risk of dying (instead of recovering).")
        st.session_state.factors_array[group][3] = c.slider("Relapsation Risk", key = f"g{group}_Relapsation Risk", min_value = 0.01, max_value = max_relapsation_factor, value = 1.3, step = 0.1, help = "Elevated relapsation rate. Higher values mean faster loss of immunity.")

_setup_subgroup_controls(st.session_state.n_subgroups)



# high_p = c1.slider("Percentage of highly susceptibles", min_value = 0.001, max_value = 0.999, value = 0.02, step = 0.001)





# c3.markdown("Highly Susceptible's \n --- ")

# infection_factor = c3.slider("Infection Risk", min_value = 0.01, max_value = max_infection_factor, value = 1.5, step = 0.1, help = "Elevated infection risk. Higher values mean higher chance of infection.")
# recovery_factor = c3.slider("Recovery Delay", min_value = 0.01, max_value = max_recovery_factor, value = 1.1, step = 0.01, help = "Delay in recovery. Higher values mean longer time to recovery.")
# death_factor = c3.slider("Death Risk", min_value = 0.01, max_value = max_death_factor, value = 1.1, step = 0.1, help = "Elevated death risk. Higher values mean increased risk of dying (instead of recovering).")
# relapsation_factor = c3.slider("Relapsation Risk", min_value = 0.01, max_value = max_relapsation_factor, value = 1.3, step = 0.1, help = "Elevated relapsation rate. Higher values mean faster loss of immunity.")


# setup general properties of the simulation

population = starting_population + starting_infectuous
endpoint = bottom_container.slider("End Timepoint", min_value = 2, max_value = 100, step = 1, value = 10)




# initialise the model
sir = model.SIR(subgroups = st.session_state.n_subgroups)

# add values to the sir model
sir.set_timespace(stop = endpoint)

sir.percentages(st.session_state.percentages)
sir.initials( values = (starting_population, starting_infectuous, 0, 0) )

sir.rates(
    infection_rate = infection_rate, 
    recovery_rate = recovery_rate,
    death_rate = theta,
    relapsation_rate = relapsation_rate
)

sir.factors(
    infection_factor = [i[0] for i in st.session_state.factors_array], 
    recovery_factor = [1/i[1] for i in st.session_state.factors_array],
    death_factor = [i[2] for i in st.session_state.factors_array],
    relapsation_factor = [i[3] for i in st.session_state.factors_array]
)

# solve the SIR Model
timespace, solutions = sir.solve()


# visualise Results

line_chart_type = line_controls.radio("Line Chart Type", ["LineChart", "PhaseChart"])

if line_chart_type == "LineChart":

    linechart = charts.LineChart(sir)
    
else:

    linechart = charts.PhaseChart(sir)

line_chart.plotly_chart(linechart, use_container_width=True)
barchart = charts.TimepointBarChart(endpoint, sir)
bar_chart.plotly_chart(barchart, use_container_width=True)