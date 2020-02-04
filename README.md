# Sparsity-Inducing Spline Penalties
This repository contains R code to implement lasso-, fusion- and mixed linear model-penalized splines. In addition, it contains code and data to generate all plots and implement the Monte Carlo simulation study for my final project in PUBH 8492 - Richly Parameterized Models, which explored the performance the three different penalization methods on various true underlying functions.  More details on the project, including a link to the resulting paper, can be found at https://sites.google.com/a/umn.edu/roland-brown/sparsity-splines.

## Contents
The contents of this repository are as follows:

  - **spline_fitting_functions.R**: This script contains all the functions necessary to estimate spline coefficients using the three competing methods, as well as functions for generating the bathtub and doppler underlying functions used in the Monte Carlo simulation study.
  - **gmst_data.csv**: Example dataset showing global mean surface temperature data.
  - **plot_gmst.R**: This script plots the Global Mean Surface Temperature Data, and shows lasso- and fusion-penalized fits along with knot coefficient values.
  - **simulation_code.R**: This script implements the Monte Carlo simulation study.
  - **Simulation Results**: This subfolder contains the simulation results, stored as .Rdata files.
  - **simulation_results.R**: This script loads results of the simulation study, and summarizes them via plots and tables.