# Introduction
The source code included is an R/Shiny application for computing inter-rater agreement statistics from simulated or real data sets. For full functionality, you'll need to install a set of tools that allow Bayesian algorithms to be run using Stan. 
If you are looking for a quicker start and are okay without the Stan computations, install the [tapModel](https://github.com/stanislavzza/tapModel) R package instead and open the app with tapModel::launchApp(). 

# Installation
Download the source code to a folder of your choice. For set up of the full functionality, you'll need:
* RTools
* Stan
* R libraries tidyverse, shiny, LaplacesDemon, and cmdstanr

# Running the app
In RStudio, open either ui.R or server.R and click the "run app" button at the top of the edit window.

# Documentation
See the [Kappa Zoo site](stanislavzza.github.io) for an explanation of how the statistics work and how to use the app. 
