#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(tidyverse)
library(LaplacesDemon) # for Mode()

source("sim_data.R")
source("load_data.R")
source("tap_analysis.R")
source("ordinal_analysis.R")
source("ta0a1p0p1_analysis.R")
source("stat_functions.R")

options(shiny.maxRequestSize = 200*1024^2) # allow large files

# Define UI for application that draws a histogram
fluidPage(
  
  tabsetPanel( id ="tabsetPanel",
               tabPanel("Load Data"),               
               tabPanel("Simulate Data"),
               tabPanel("t-a-p"),
               tabPanel("Ordinal t-a-p"),
               tabPanel("t-a0,a1-p0,p1")
               
  ),
  
  conditionalPanel(
    condition = "input.tabsetPanel == 'Load Data'",
    data_load_ui("load_data")
  ),
  
  conditionalPanel(
    condition = "input.tabsetPanel == 'Simulate Data'",
    data_sim_ui("sim_data")
  ),
  
  conditionalPanel(
    condition = "input.tabsetPanel == 't-a-p'",
    tap_analysis_ui("t-a-p")
  ),
  
  conditionalPanel(
    condition = "input.tabsetPanel == 'Ordinal t-a-p'",
    ordinal_analysis_ui("ordinal")
  ),
  
  conditionalPanel(
    condition = "input.tabsetPanel == 't-a0,a1-p0,p1'",
    ta0a1p0p1_analysis_ui("t-a0a1-p0p1")
  )
  
  
)
