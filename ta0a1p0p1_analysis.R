ta0a1p0p1_analysis_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel( 
      uiOutput(ns("targetVar")), 
      uiOutput(ns("targetRangeIn")),
      uiOutput(ns("action_tap")),
      checkboxInput(ns("p0p1"),"Include p0,p1")
    ),
    
    mainPanel(
      h1("t-a0,a1-p0,p1 parameter estimates"),
      plotOutput(ns("mcmc_distro"),height = 350),
      
    )
  )
}

ta0a1p0p1_analysis_server <- function(id, ratings) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    rating_counts <- reactiveVal(NULL)
    param_hat     <- reactiveVal(NULL)
    mcmc_results  <- reactiveVal(NULL)
    mcmc          <- reactiveVal(NULL)
    mcmc_distro   <- reactiveVal(NULL)
    
    #INPUT--------------------- Select the Target Var to predict-----------------------
    output$targetVar <- renderUI({ 
      # turning off the newvars update for now--it's quite annoying
      #rvals$newVarExists # update the list of names if a new var is created
      req(ratings())
      namelist <- as.list(sort(names(ratings())))
      namelist <- namelist[namelist != "SubjectID__"]
      select <- isolate(input$targetVar) # save the old choice if there is one
      selectInput(ns("targetVar"), "Ratings",  choices = namelist , selected = select )
    })
    
    #INPUT -------------------- Select the target inclass range--------------------------------
    output$targetRangeIn <- renderUI({ 
      req(ratings(), input$targetVar)
      
      namelist <- as.list(names(table(ratings()[[input$targetVar]])))
      
      selectInput(ns("targetRangeIn"), "Choose the inclass range", namelist, multiple = TRUE )
    })
    
    #INPUT -------------------- Button to compute the tap stats
    output$action_tap <- renderUI({
      req(input$targetRangeIn)
      
      actionButton(ns("action_tap"),"Compute")
    })
    
    observeEvent(input$action_tap, {
      
      # ratings are in df, but need to be filtered
      tdf <-  ratings() |> 
        select(SubjectID__, rating = input$targetVar) |> 
        na.omit() |> 
        group_by(SubjectID__) |>
        summarise(N_r = n(),
                  N_c = sum(rating %in% input$targetRangeIn)) |> 
        filter(N_r > 1) # need at least pairs of raters
      
      # set the reactive values to the output   
      rating_counts(tdf) # set the counts
      param_hat( iterative_optim(rating_counts()) ) # param estimates
      
      mcmc_results(  get_mcmc_results2(rating_counts(), input$p0p1) ) 
      mcmc(mcmc_results()$stats)
      mcmc_distro(mcmc_results()$draws)
    })
    
    output$mcmc_distro <- renderPlot({
      
      req(mcmc_distro())
      
      #input$action_tap2 # react to it
      plot_draw_densities(mcmc_distro() )
      
    })
  })
}