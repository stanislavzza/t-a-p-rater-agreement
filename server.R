#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define server logic required to draw a histogram
function(input, output, session) {
  
  ratings <- reactiveVal()
  
  # ratings from either a csv or simulation
  csv_data <- data_load_server("load_data")
  sim_data <- data_sim_server("sim_data")
  
  # update ratings from either source
  observeEvent(csv_data(), {
    ratings(csv_data())
  })
  
  observeEvent(sim_data(), {
    ratings(sim_data())
  })
  
  # run analysis
  tap_analysis_server("t-a-p", ratings)
  ordinal_analysis_server("ordinal", ratings)
  ta0a1p0p1_analysis_server("t-a0a1-p0p1", ratings)

}
