# Effect of population heterogeneity in disease dynamics
### A project for Computational Epidemiology (Course: 467294-HS2021)

[![Open the interactive Web-App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/noahhenrikkleinschmidt/computational_epidemiology_sir_model/main/main.py)

This project models the effect of two subgroups within a population. These subgroups are colloquially termed either "normally" susceptible ($S_n$, default average person) or "highly" susceptible ($S_h$, a subgroup that differs from the average by some parameters). $S_h$ are modelled as a linear-transform of $S_n$ and may differ in any combination of infection rate, recovery rate, death rate, and/or relapsation rate (loss of immunity) from $S_n$ by a scalar increase or decrease relative to $S_n$'s respective transition rate.


## The SIR Model
---

### Susceptible Population $S$ 

This project assumes the total population $S$ consists of two groups as per:
	

$$
S = S_n + S_h = \alpha_n  \cdot S + \alpha_h \cdot S
$$

where $\alpha_i$ are the percentages of each sub-group within the population, $\alpha_n + \alpha_h = 1$. 

### Infectious Population $I$ 

As both $S_n$ and $S_h$ are assumed to have different susceptibilities, they are assigned two different rates of infection: $\beta_n$ and $\beta_h$, where we assume for simplicity sake (see Introduction) that $\beta_h = k\cdot\beta_n$ for some $k \in \mathbb{R}$. Given the scalar relation between $\beta_n$ and $\beta_h$ we abbreviate $\beta_n \equiv \beta$, as we will do for any transition rate henceforth.


### Removed Population $R$ 

The removed population $R$ consists of two sub-populations, namely the recovered population $R'$ and the deceased population $D$. Infectious individuals may transition into either $R'$ or $D$, given two separate transition rates. 

We define the recovery rate $\gamma$ and death rate $\theta$ whose scaled rates for the "highly" susceptibles are $\gamma_h = j\cdot\gamma$ and $\theta_h = q \cdot \theta$ for some $j, q \in \mathbb{R}$, respectively.

We further define a relapsation rate $\delta$ that mimics how immunity is lost over time (perhaps due to newly emerging mutant pathogens) and allows individuals from $R'$ to relapse into $S$. Here as well "highly" susceptibles are linearly scaled as $\delta_h = h \cdot \delta$ for some $h \in \mathbb{R}$. 

### System of Equations 

The outlined model allows transitions between Susceptible $S$ to Infectious $I$, from $I$ to either Removed population $R'$ or $D$, and from $R'$ back to $S$. This allows the following set of first order differential equations to be formulated:

$$
\begin{align}
\frac{dS(t)}{dt} = {}&  -\beta \Phi_k \cdot S(t) \cdot I(t) + \delta\Phi_h \cdot R'(t) \\

\frac{dI(t)}{dt} = {}&  \beta \Phi_k \cdot S(t) \cdot I(t) - \gamma\Phi_j \cdot I(t) - \theta\Phi_q\cdot I(t) \\

\frac{dR'(t)}{dt} = {}& \gamma\Phi_j \cdot I(t) - \delta\Phi_h \cdot R'(t) \\

\frac{dD(t)}{dt} = {}&  \theta\Phi_q\cdot I(t)
\end{align}
$$

Where $\Phi_x = \alpha_n + x\cdot \alpha_h$ describes the impact of "highly" susceptibles within the total population for some scalar factor $x$ of a given transition. 


This set of differential equations now has four separate rates of transition: (1) $\beta$, the rate of infection, where "highly" susceptibles are distinguished by a $k$-fold increase/ decrease. (2) $\gamma$, the rate of recovery, where "highly" susceptibles are distinguished by a $j$-fold increase/decrease. (3) $\delta$, the relapsation rate, where "highly" susceptibles are distinguished by an $h$-fold increase/drecrease. And (4) $\theta$, the death rate, where "highly" susceptibles are distinguished by a $q$-fold increase/decrease compared to "normally" susceptibles. All of these are editable to explore different scenarios of heterogeneous populations. 
