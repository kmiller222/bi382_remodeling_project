# bi382_remodeling_project
### SIRV disease model of the impacts of vaccination on a hypothetical influenza outbreak, based on a recreation of the paper, Impact of influenza vaccine-modified infectivity on attack rate, case fatality ratio, and mortality (Nah et al. 2020)

This project is my recreation of the modeling of influenza in a 2020 paper titled “Impact of influenza vaccine-modified infectivity on attack rate, case fatality ratio and mortality” by Nah et al., which I completed for my Ecological Modeling (BI382) class at Colby College in 2020.

## Model Code
The following code is what I used as a base SIRV model. The SIRV model uses differential equations in order to relate continuous rates of change of these population classes to the biological processes of disease transmission, which I created in R by using the `deSolve` package. The SIRV model adds classes of susceptible vaccinated (V) and infected vaccinated (Iv), in addition to the usual susceptible (S), infected (I), and recovered (R) population groups of an SIR model to also incorporate vaccine-modified effects. Parameters were changed throughout the project to reflect different vaccination and disease-characteristic scenarios. 

```
library(deSolve)

# set parameters
N0 <- 1000
R0 <- 1.5 # basic reproductive number
gamma <- 0.33  # recovery rate of non-vaccinated infecteds per infected per day
delta <- 0.0 # disease-induced death for non-vaccinated individuals per infected per day
ms <- 0.6  # modification of susceptibility for vaccinated individuals
md <- 0.5  # modification of death rate for vaccinated individuals
steps <- 100
mi <- 0.5
mr <- 0.5
rho <- 0.5

# SIRV model
SIRV <- function(t, init, parms) {
  with(as.list(c(init, parms)), {
    dS <- -beta*S*(I + mi*Iv)
    dV <- -beta*ms*V*(I + mi*Iv)
    dI <- beta*S*(I + mi*Iv) - gamma*I - delta*I
    dIv <- beta*ms*V*(I + mi*Iv) - mr*gamma*Iv - md*delta*Iv
    dR <- gamma*I + mr*gamma*Iv 
    return(list(c(dS, dV, dI, dIv, dR)))
  })
}

initial_infections <- 20
S0 <- (1-rho)*N0 - initial_infections*(1-rho)
V0 <- rho*N0 - initial_infections*rho
model_days <- 100
times <- seq(from = 1, to = model_days, by = 0.25)

sirv_parms <- c(beta = (R0*((delta + gamma)/N0)), mi=mi, ms=ms, mr=mr, md=md, 
                    gamma=gamma, delta=delta)
sirv_init <- c(S=S0, V=V0, I=initial_infections*(1-rho), 
               Iv=initial_infections*rho, R=0)

sirv_out <- ode(func = SIRV, parms = sirv_parms, times = times, y = sirv_init)
sirv_result = as.data.frame(sirv_out)
```


## References

Nah, K., Alavinejad, M., Rahman, A., Heffernan, J. M., & Wu, J. (2020). Impact of influenza vaccine-modified infectivity on attack rate, case fatality ratio and mortality. Journal of Theoretical Biology, 492, 110190.

Moore, C. (2020). Disease Modeling in R [BI382 Colby College, Class Notes].
