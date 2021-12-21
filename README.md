# Effect of population heterogeneity in disease dynamics
### A project for Computational Epidemiology (Course: 467294-HS2021)

This project models the effect of two subgroups within a population that are either "normally" susceptible (default average person) or "highly" susceptible to a disease (perhaps due to immuno-suppression because of other medical issues). As such, "normally" susceptibles ($S_n$) and "highly" susceptibles ($S_h$) are assumed to have also different rates of recovery.

## SIR Model 
---
### Susceptible Population $S$
This project assumes the total population $S$ consists of two groups as per:

$$
S = S_n + S_h = \alpha_n  \cdot S + \alpha_h \cdot S = \alpha_n  \cdot S + (1-\alpha_n)\cdot S
$$

where $\alpha_i$ are the percentages of each group within the population. 

### Infectious Population $I$
As both $S_n$ and $S_h$ are assumed to have different susceptibilities, they are assigned two different rates of infection: $\beta_n$ and $\beta_h$, where we assume for simplicity sake that $\beta_h = k\cdot\beta_n$ for some $k > 1 \in \R$. 


### Removed Population $R$
Removed (recovered or deceased) cannot be infected nor infect others. Both $S_n$ and $S_h$ have a recovery rate that assigns them to this population. Yet, we assume that the two groups have different rates of recovery $\gamma_n$ and $\gamma_h$. We again assume for simplicity sake that $\gamma_h = j\cdot\gamma_n$ for some $j < 1 \in \R$. The total rate at which $R$ increases is hence $\gamma_n + \gamma_h = \gamma_n + j\gamma_n = (1+j)\gamma_n$.

### Preliminary SIR System of Equations
As both infection rate and recovery rate are ultimately just determined by the rate for the "normally" susceptibles we abbreviate $\beta_n \equiv \beta$ and $\gamma_n \equiv \gamma$.

The implemented system of differential equations is hence: 

$$
\begin{align}

\frac{dS(t)}{dt} = {}&  [\ -\beta \cdot S_n(t) - k\beta \cdot S_h(t) \ ] \cdot I(t) \\
\frac{dI(t)}{dt} = {}&  [\ \beta \cdot S_n(t) + k\beta \cdot S_h(t) \ ] \cdot I(t) - (1 + j)\gamma\cdot I(t) \\
\frac{dR(t)}{dt} = {}&  (1 + j)\gamma\cdot I(t)

\end{align}
$$


This system we can simplify slightly by aggregating the infectous transition for $S_n$ and $S_h$ as follows, thus modifying $\frac{dS(t)}{dt}$ and $\frac{dI(t)}{dt}$:

$$
\begin{align}

\frac{dS(t)}{dt} = {}&  -\beta(\alpha_n  + k (1 -\alpha_n)) \cdot S(t)\cdot I(t) \\
\frac{dI(t)}{dt} = {}&  \beta(\alpha_n  + k (1 -\alpha_n)) \cdot S(t)\cdot I(t) - (1 + j)\gamma\cdot I(t) \\

\end{align}
$$

### Relapsation of susceptibles
In our preliminary model we assume two groups within the population, the "normal" and "highly" susceptibles. We further assume that, although they both have different rates of infection and recovery, they ultimately endure the same course of disease, going from susceptible, to infectous, and ultimately to removed. However, what if this assumption were not true, and a small part of the population were to relapse into the susceptible population $S$ after having endured the infectuous stage? 

This could be the case for immuno-suppressed inviduals who did survive the disease due to medical care but were unable to fully develop their own immune response. Alternatively, normally susceptibles may have endured a more severe course of disease that weakened their immune system, thus pushing them to the highly susceptibles after having left the infectuous stage. 

#### Relapsation Rate $\delta$
Both of the above scenarious are somewhat complex, and shall be addressed in a simplified manner by introducing a relapsation rate $\delta_i$ for both normally and highly susceptibles. $\delta_i$ allows removed people to leave $R$ to re-enter $S$ as part of the highly susceptibles $S_h$. Thereby, $\frac{dS(t)}{dt}$ will be modified to receive an additive term from $R$, and $\frac{dR(t)}{dt}$ will receive a corresponding subtractive term. As we assume normal and highly susceptibles will have different rates of relapsation, we will use our default asumption of $\delta_h = h \cdot \delta_n$ for some $h \in \R$ (no restriction to bigger or smaller 1).

### Relapsation SIR System of Equations
We retain $\beta_n \equiv \beta$ and $\gamma_n \equiv \gamma$, and abbreviate $\delta_n \equiv \delta$.
The implemented system of differential equations is hence: 

$$
\begin{align}

\frac{dS(t)}{dt} = {}&  -\beta(\alpha_n  + k (1 -\alpha_n)) \cdot S(t) \cdot I(t) + \delta(1 + h) \cdot R(t) \\
\frac{dI(t)}{dt} = {}&  \beta(\alpha_n  + k (1 -\alpha_n)) \cdot S(t) \cdot I(t) - (1 + j)\gamma\cdot I(t) \\

\frac{dR(t)}{dt} = {}&  (1 + j)\gamma\cdot I(t) - \delta(1 + h) \cdot R(t)

\end{align}
$$


### Lethality versus Recovery
The model outlined so far has treated recovered and deceased as equivalent. However, considering the assumtion above of some more susceptible groups that may also have a higher relapsation rate, the concept of lethality should be taken into consideration as well. Immuno-suppressed individuals may very likely have a higher risk of actually dying rather than recovering. This would impact the subtractive term of $\frac{dI(t)}{dt}$ (additive in $\frac{dR(t)}{dt}$), as only a sub-group of removed may be able to relapse, namely the ones who survived. To address this we modify the model to split the current $R$ into two groups $R'$ (recovered, who may relapse) and $D$ (deceased, who may not relapse). Admitedly, this does leave the _SIR_ model for an "_SIRD_" model (not sure about terminlolgy?), but seems to mimic disease dynamics better under the given assumptive parameters. 

#### Death rate $\theta$
To split the current $R$ (removed) we adjust the rate of transition from $I$ to $R$ for the two new subgroups of $R$, we introduce a lethality rate $\theta$ as a percentage of people who transition to $R$ due to death rather than recovery. $\theta$ modifies the current recovery rate $\gamma$ to a new recovery rate $\gamma' = (1 - \theta)\gamma$. 
We again assume that highly susceptibles have a different death rate than normally susceptibles: $\theta_h = q \cdot \theta_n$ for some $q > 1 \in \R$. 


The new transition from $I$ to $R'$ or $D$ is hence: 

$$
\begin{align}
\frac{dI(t)}{dt} = {}&  \beta(\alpha_n  + k (1 -\alpha_n)) \cdot S(t) \cdot I(t) - (1 + j)\gamma'\cdot I(t) - (1 + q)\theta\cdot I(t) \\

\frac{dR'(t)}{dt} = {}&  (1 + j)\gamma'\cdot I(t) - \delta(1 + h) \cdot R'(t) \\
\frac{dD(t)}{dt} = {}&  (1 + q)\theta\cdot I(t)

\end{align}
$$


## The final SIR(D) Model
The final implementation allows for heterogeneity within the population on infection level, recovery level, as well as relapsation level. Two groups within the population are considered, "normally susceptible" and "highly susceptible" people. Highly susceptibles are destinguished as a linearly transformed version of normal susceptibles. To model lethality impacts and relapsation rates of the different subgroups, the classical $R$ (removed) group was split into a recovered group ($R'$) and a deceased group $D$. 

The final system of differential equations hence looks as follows: 

$$
\begin{align}

\frac{dS(t)}{dt} = {}&  -\beta(\alpha_n  + k (1 -\alpha_n)) \cdot S(t) \cdot I(t) + \delta(1 + h) \cdot R'(t) \\
\frac{dI(t)}{dt} = {}&  \beta(\alpha_n  + k (1 -\alpha_n)) \cdot S(t) \cdot I(t) - (1 + j)\gamma'\cdot I(t) - (1 + q)\theta\cdot I(t) \\

\frac{dR'(t)}{dt} = {}&  (1 + j)\gamma'\cdot I(t) - \delta(1 + h) \cdot R'(t) \\
\frac{dD(t)}{dt} = {}&  (1 + q)\theta\cdot I(t)

\end{align}
$$

This set of differential equations now has four separate rates of transition: (1) $\beta$, the rate of infection, where highly susceptibles are distinguished by a $k$-fold increase. (2) $\gamma$, the rate of recovery, where highly susceptibles are distinguished by a $j$-fold decrease. (3) $\delta$, the relapsation rate, where highly susceptibles are distinguished by a $h$-fold change compared to normally susceptibles. And (4) $\theta$, the death rate, where highly susceptibles are distinguished by a $q$-fold increase compared to normally susceptibles.
¨

> Alright, so far so good, let's try to implement this into code >_< ...
