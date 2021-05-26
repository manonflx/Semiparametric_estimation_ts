# Semidiscrete_estimation_ts
R code to reproduce all plots/tables in the master's thesis: "Semiparametric estimation for time series: a frequency domain approach based on optimal transportation theory"

Two functions are mainly used: 

1. Function whittle $=$ function(theta, phi, $\mathrm{H}, \mathrm{p}, \mathrm{q}, \mathrm{y}$ ) where:
- theta is a vector of length $p$.
- $p h i$ is a vector of length $q$.
- $H$ self-similarity parameter (recap: ARIMA(0,d,0) process $d=H-1 / 2$ )
- $p, q$ are optional integers specifying the $A R$ and $M A$ orders of the FARIMA model
- $y$ is a numeric vector representing a time series
