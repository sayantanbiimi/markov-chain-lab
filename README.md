# Markov Chain Lab

An interactive Shiny app for teaching recurrent and transient states.

Copyright: Sayantan Banerjee, IIM Indore

## Run locally

```r
install.packages("shiny")
shiny::runApp("markov-chain-lab")
```

The app uses only the `shiny` package. It includes four finite-state examples,
exact structural classification, sample paths, finite-horizon return
probabilities, transition-matrix powers, and a simulation of the biased random
walk on the integers.

## Classroom activities

1. In the “one transient state” example, keep rerunning paths from state 1.
   Can a path return even though the state is transient?
2. Compare the periodic chain with the other recurrent chain. At which times
   can it return?
3. Set the random-walk probability to 0.50 and then 0.55. Why does a tiny drift
   change the classification even though short paths look similar?
4. Increase the simulation horizon. Explain why the simulated return fraction
   estimates a finite-horizon probability, not automatically the eventual one.
