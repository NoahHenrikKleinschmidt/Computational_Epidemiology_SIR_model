# Effect of population heterogeneity in disease dynamics
### A project for Computational Epidemiology (Course: 467294-HS2021)

[![Open the interactive Web-App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/noahhenrikkleinschmidt/computational_epidemiology_sir_model/main/main.py)

This project models the effect of a variety of subgroups within a population on disease dynamics. To that end an SIR model was developed that models a susceptible population of various subgroups that differ in their infection rate, recovery rate, death rate, and/or relapsation rate (loss of immunity). Four different populations are assumed to that end: Susceptible $S$, Infectuous $I$, Recovered $R$, and Deceased $D$. Colloqually the term "population" will be used to refer to any or all non-deceased populations.

## The SIRD Model
---

This model assumes the population contains a reference subgroup termed $S_n$ (for "normally" susceptibles) that transitions between Susceptible population $S$, to Infectuous $I$, Recovered $R$, or Desceased $D$ with some default transition rates (see below). This model allows for an arbitrary number of divergent subgroups $S_i$ within the population $S$, represented by a percentage $\alpha_i$. 

Consequently, it holds that $\alpha_n + \sum \alpha_i = 1$, where $\alpha_n$ is the percentage of the reference subgroup within $S$. And the total population is thus given by:

$$
S = S_n + \sum S_i = \alpha_n \cdot S + \sum (\alpha_i  \cdot S) = (\alpha_n + \sum \alpha_i)\cdot S
$$

### Divergent Subgroups

For the purpose of this model each divergent subgroup is modelled as a linear-transform of the reference group. Hence, each $S_i$ is distinguished by a scalar factor $x$ for each transition rate. As this model considers four transition rates (see below) each $S_i$ will be associated with a tuple of four corresponding $x_t \in \mathbb{R}$, one for each transition $t$. The impact of each divergent subgroup is given by its percentual representation within the overall population. Consequently, we can define a weight function $\Phi_t$ that assignes the impact of each divergent subgroup for a given transition $t$

$$
\Phi_t = \alpha_n + \sum_{i}(x_{i,t}\cdot \alpha_i)
$$

where $x_{i,t}$ is the respective scalar factor for a given transition $t$ of subgroup $i$. 
The actual transition rate for the overall system $\lambda_t$ is hence given by

$$
\lambda_t = \Phi_t \cdot t
$$

where $t \in \mathbb{R}$ denotes the default transition rate. 


### Transitions
This model allows for four transitions between the different populations. 
(1) Susceptibles may transition into Infectuous with a default infection rate $\beta$. 
(2) Infectuous may transition into Recovered with a default recovery rate $\gamma$. 
(3) Alternatively, Infectuous may transition into Deceased with a default death rate $\theta$. 
(4) Recovered may transition back into Susceptibles when they loose their immunity over time with a default relapsation rate $\delta$. 


### System of Equations 

The outlined model allows the following set of first order differential equations to be formulated:

$$
\begin{align}
\frac{dS(t)}{dt} = {}&  -\lambda_\beta \cdot S(t) \cdot I(t) +\lambda_\delta \cdot R(t) \\

\frac{dI(t)}{dt} = {}&  \lambda_\beta \cdot S(t) \cdot I(t) - (\lambda_\gamma + \lambda_\theta )\cdot I(t) \\

\frac{dR(t)}{dt} = {}& \lambda_\gamma \cdot I(t) - \lambda_\delta \cdot R(t) \\

\frac{dD(t)}{dt} = {}&  \lambda_\theta \cdot I(t)
\end{align}
$$


# Old stuff from here on...


Where $\Phi_x = \alpha_n + x\cdot \alpha_h + x'\cdot \alpha_l$ describes the impact of "highly" susceptibles within the total population for some scalar factors $x$ and $x'$ of a given transition for the two non-average subgroups. 


This set of differential equations now has four separate rates of transition: (1) $\beta$, the rate of infection, where "highly" susceptibles are distinguished by a $k$-fold increase/ decrease. (2) $\gamma$, the rate of recovery, where "highly" susceptibles are distinguished by a $j$-fold increase/decrease. (3) $\delta$, the relapsation rate, where "highly" susceptibles are distinguished by an $h$-fold increase/drecrease. And (4) $\theta$, the death rate, where "highly" susceptibles are distinguished by a $q$-fold increase/decrease compared to "normally" susceptibles. All of these are editable to explore different scenarios of heterogeneous populations. 





These subgroups are colloquially termed either "normally" susceptible ($S_n$, default average person) or "highly" susceptible ($S_h$, a subgroup that differs from the average by some parameters such as an elevated infection rate) and "lowly" susceptible ($S_l$, that differs from the average by some parameters such as elevated recovery rate). $S_h$ and $S_l$ are modelled as a linear-transform of $S_n$ and may differ in any combination of infection rate, recovery rate, death rate, and/or relapsation rate (loss of immunity) from $S_n$ by a scalar increase or decrease relative to $S_n$'s respective transition rate. 

### Susceptible Population $S$ 

This project assumes the total population $S$ consists of two groups as per:
	

$$
S = S_n + S_h = \alpha_n  \cdot S + \alpha_h \cdot S
$$

where $\alpha_i$ are the percentages of each sub-group within the population, $\alpha_n + \alpha_h = 1$. 

### Infectious Population $I$ 

As all $S_n$, $S_h$, and $S_l$ are assumed to have different susceptibilities, they are assigned three different rates of infection: $\beta_n$, $\beta_h$,  and $\beta_l$, where we assume for simplicity sake (see Introduction) that $\beta_h = k\cdot\beta_n$ and $\beta_l = k'\cdot\beta_n$ for some $k, k' \in \mathbb{R}$. Given the scalar relation between $\beta_n$, $\beta_h$, and $\beta_l$ we abbreviate $\beta_n \equiv \beta$, as we will do for any transition rate henceforth.


### Removed Population $R$ 

The removed population $R$ consists of two sub-populations, namely the recovered population $R'$ and the deceased population $D$. Infectious individuals may transition into either $R'$ or $D$, given two separate transition rates. 

We define the recovery rate $\gamma$ and death rate $\theta$ whose scaled rates for the "highly" susceptibles are $\gamma_h = j\cdot\gamma$ and $\theta_h = q \cdot \theta$ and $\gamma_l = j'\cdot\gamma$ and $\theta_l = q' \cdot \theta$ for the "lowly" susceptibles, for some $j, q, j', q' \in \mathbb{R}$, respectively.

We further define a relapsation rate $\delta$ that mimics how immunity is lost over time (perhaps due to newly emerging mutant pathogens) and allows individuals from $R'$ to relapse into $S$. Here as well "highly" susceptibles are linearly scaled as $\delta_h = h \cdot \delta$ for some $h \in \mathbb{R}$, as are the "lowly" susceptibles by $\delta_l = h' \cdot \delta$ for some $h' \in \mathbb{R}$. 
