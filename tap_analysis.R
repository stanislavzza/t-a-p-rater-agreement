tap_analysis_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel( 
      uiOutput(ns("targetVar")), 
      uiOutput(ns("targetRangeIn")),
      uiOutput(ns("targetRangeOut")),
      uiOutput(ns("distro")), 
      uiOutput(ns("N_r")),
      uiOutput(ns("scale")),
      uiOutput(ns("xaxis")),
      uiOutput(ns("mcmc")),
      uiOutput(ns("action_tap"))
    ),
    
    mainPanel(
      h1("t-a-p parameter estimates"),
      plotOutput(ns("tap"),height = 350),
      plotOutput(ns("ll_path"),height = 350),
      plotOutput(ns("mcmc_distro"),height = 350),
      
    )
  )
}

tap_analysis_server <- function(id, ratings) {
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
    
    #INPUT -------------------- Select the target outclass range--------------------------------
    output$targetRangeOut <- renderUI({ 
      req(ratings(), input$targetRangeIn)
      
      namelist <- as.list(names(table(ratings()[[input$targetVar]])))
      # remove the ones in input$targetRangeIn
      namelist <- namelist[!namelist %in% input$targetRangeIn]
      
      selectInput(ns("targetRangeOut"), "Choose the outclass range", namelist, multiple = TRUE, selected = namelist )
    })
    
    #INPUT -------------------- Button to compute the tap stats
    output$action_tap <- renderUI({
      req(ratings(), input$targetVar, input$targetRangeIn,input$targetRangeOut)
      if (length(input$targetRangeIn) == 0 ) return (NULL)
      if (length(input$targetRangeOut) == 0 ) return (NULL)
      
      actionButton(ns("action_tap"),"Compute")
    })
    
    #INPUT -------------------- t-a-p display options: combined distribution
    output$N_r <- renderUI({
      req(rating_counts()) # computed after the button is pressed
      
      choice_list <- rating_counts() |> 
        count(N_r) |> 
        slice_max(n = 10, order_by = n) |> 
        mutate(label = str_c(N_r," (n = ",n,")")) |> 
        select(label, N_r) |> 
        deframe() |> 
        as.list()
      
      selectInput(inputId = ns("N_r"),
                  label = "Number of raters",
                  choices = choice_list,
                  multiple = FALSE, 
                  selectize = FALSE)
    })
    
    output$distro <- renderUI({
      req(ratings())
      
      checkboxInput(ns("distro"),
                    "Show combined distribution",
                    value = TRUE)
    })
    
    #INPUT -------------------- t-a-p display options: rescale by t 
    output$scale <- renderUI({
      req(ratings())
      
      checkboxInput(ns("scale"),
                    "Scale density by t",
                    value = TRUE)
    })
    
    #INPUT -------------------- t-a-p display options: compute MCMC?
    output$mcmc <- renderUI({
      req(ratings())
      
      checkboxInput(ns("mcmc"),
                    "Compute MCMC (slow)",
                    value = FALSE)
    })
    
    #INPUT -------------------- t-a-p display options: choose x-axis
    output$xaxis<- renderUI({
      req(ratings())
      
      radioButtons(ns("xaxis"),
                   "LL x-axis",
                   c("t" = "t", "a" = "a", "p" = "p"),
                   selected = "t",
                   inline = TRUE)
    })
    
    # PROCESS DATA -----------------------
    observeEvent(input$action_tap, {
      
      target_ranges <- c(input$targetRangeIn,input$targetRangeOut)
      
      # ratings are in df, but need to be filtered
      tdf <-  ratings() |> 
        select(SubjectID__, rating = input$targetVar) |> 
        na.omit() |> 
        filter(rating %in% target_ranges) |>
        group_by(SubjectID__) |>
        summarise(N_r = n(),
                  N_c = sum(rating %in% input$targetRangeIn)) |> 
        filter(N_r > 1) # need at least pairs of raters
      
      # set the reactive values to the output   
      rating_counts(tdf) # set the counts
      param_hat( iterative_optim(rating_counts()) ) # param estimates
    
      # do mcmc if box is checked
      if(isolate(input$mcmc)) {
        mcmc_results(  get_mcmc_results(rating_counts()) ) # reactive val
        mcmc(mcmc_results()$stats)
        mcmc_distro(mcmc_results()$draws)
      }

    })
    
    
    # plot distributions of x and y
    output$tap<- renderPlot({
      
      req(param_hat())
      
      input$action_tap # refresh when action
      
      # estimated parameters (hat)
      t_hat <- param_hat()$t
      a_hat <- param_hat()$a
      p_hat <- param_hat()$p
      
      # find the mode for raters per subject
      N_r <- isolate(input$N_r)
      
      # if it hasn't been set yet, set it to the max n value
      if(is.null(N_r)) {
        N_r <- rating_counts() |> 
          count(N_r) |> 
          arrange(desc(n)) |> 
          slice(1) |> 
          pull(N_r) 
      } else {
        N_r <- as.integer(N_r)
      }
      
      counts <- rating_counts()|> 
        filter(N_r == !!N_r) 
      
      # generate the title with the parameters
      my_title <- str_c("Estimated t = ", round(t_hat, 2), 
                        ", a = ", round(a_hat, 2), 
                        ", p = ", round(p_hat, 2))
      
      # model distribution parameters
      T1_mu <- a_hat + (1-a_hat)*p_hat
      T0_mu <- (1-a_hat)*p_hat
      
      
      ################# Switch presentation depending on checkbox
      if(!input$distro){
        
        # modeled probability values for each case
        pdf_0 <- tibble(T_i = 0, N_c = 0:N_r,
                        prob = dbinom(N_c, N_r, T0_mu),
                        Distribution = "Model")
        
        pdf_1 <- tibble(T_i = 1, N_c = 0:N_r,
                        prob = dbinom(N_c, N_r, T1_mu),
                        Distribution = "Model")
        
        pdf <- rbind(pdf_0, pdf_1) 
        
        if(input$scale == TRUE){
          pdf <- pdf |> 
            mutate(prob = if_else(T_i == 1, prob*t_hat, prob*(1-t_hat)))
          
        }
        # make the scale nicer 
        pdf <- pdf |> 
          mutate(T_i = factor(T_i, labels = c("Out-class","In-class")))
        
        # plot the probability of C = 1 for both cases
        pdf |> 
          ggplot(aes(x = N_c, y = prob, group = T_i, 
                     color = T_i)) +
          geom_line(linewidth = 1, linetype = "dashed") +
          # set colors manually
          scale_color_manual(values = c("black", "#096DB9")) +
          theme_bw() +
          ggtitle(my_title) +
          theme(text=element_text(size=20))
      } else {
        
        # empirical rating distribution
        true_dist <- counts %>%
          count(N_c) %>%
          mutate(prob = n/sum(n)) |> 
          select(-n)
        
        # model values
        pdf <- tibble(N_c = 0:N_r,
                      prob = t_hat*dbinom(N_c, N_r, T1_mu) +
                        (1-t_hat)*dbinom(N_c, N_r, T0_mu))
        
        
        pdf |> 
          ggplot(aes(x = N_c, y = prob)) +
          geom_line(linewidth = 1, color = "#096DB9", linetype = "dashed") +
          # draw line up from x-axis
          geom_segment(data = true_dist, aes(x = N_c, xend = N_c, y = 0, yend = prob),
                       linetype = "solid", color = "#222222") +
          geom_point(data = true_dist, aes(x = N_c, y = prob),
                     color = "#222222", size = 2) +
          theme_bw() +
          ggtitle(my_title)+
          theme(text=element_text(size=15)) +
          ylab("Probability") +
          xlab("Rating agreement for in-class") +
          coord_cartesian(clip = 'off')
      }
      
    })
    
    # plot the log-likelihood path
    output$ll_path <- renderPlot({
      
      # if(is.null(rvals$counts)) return(NULL)
      req(rating_counts())
    
      input$action_tap # react to it
      
      # set up the plotting frame on [0,1]
      lls <- tibble(x = seq(0,1,.01))
      
      # a hack to overcome a problem with geom_ribbon, which leaves gaps
      # between regions. Need to add more points at the breaks 
      if(input$mcmc & !is.null(mcmc())) {
        shading <- mcmc() |>
          filter(var == input$xaxis) |> 
          as.list()
        
        spread_point <- function(x) seq(x - .008, x + .008, .002)
        
        new_x <- data.frame(x = c(spread_point(shading$p05),
                                  spread_point(shading$p25),
                                  spread_point(shading$median-.01),
                                  spread_point(shading$median+.01),
                                  spread_point(shading$p75), 
                                  spread_point(shading$p95)))
        lls <- bind_rows(lls, new_x) |> 
          arrange(x)
        
      }
      
      # get the tap parmeters
      t_hat <- param_hat()$t
      a_hat <- param_hat()$a
      p_hat <- param_hat()$p
      ll_hat <- param_hat()$ll
      
      # fleiss kappa
      kappa_stats <- fleiss_kappa(rating_counts())
      ll_kappa <- kappa_stats$ll
      
      if(input$xaxis == "t") {
        lls <- lls |> 
          rowwise() |>
          mutate(ll = -log_likelihood(c(x,a_hat,p_hat), rating_counts()))
        
        kappa_intercept <- kappa_stats$t
        model_intercept <- t_hat
      }
      
      if(input$xaxis == "a") {
        lls <- lls |> 
          rowwise() |>
          mutate(ll = -log_likelihood(c(t_hat,x,p_hat), rating_counts()))
        
        kappa_intercept <- kappa_stats$a
        model_intercept <- a_hat
      }
      
      if(input$xaxis == "p") {
        lls <- lls |> 
          rowwise() |>
          mutate(ll = -log_likelihood(c(t_hat,a_hat,x), rating_counts()))
        
        kappa_intercept <- kappa_stats$p
        model_intercept <- p_hat
      }
      
      # plotting, starting with mcmc if available.
      if(isolate(isTruthy(input$mcmc))) { #  & nrow(mcmc_distro()) > 0
        
        lls <- lls |> 
          mutate(region = case_when(x < shading$p05 ~ 0,
                                    x < shading$p25 ~ 1,
                                    x < shading$median - .01 ~ 2,
                                    x < shading$median + .01 ~ 3,
                                    x < shading$p75 ~ 4,
                                    x < shading$p95 ~ 5, 
                                    TRUE ~ 6),
                 region = as.factor(region))
        
        # start the plot    
        gg <- ggplot(lls, aes(x = x, y = ll))  +
          geom_ribbon( aes( fill = region, ymax = ll), 
                       ymin = -Inf, 
                       position = position_identity()) +
          scale_fill_manual(values = c("#FFFFFF00", 
                                       "#BBBBFF55", 
                                       "#4444FF88", 
                                       "#1111FF",
                                       "#4444FF88", 
                                       "#BBBBFF55", 
                                       "#FFFFFF00"),
                            drop = TRUE,  # omit unused factors
                            limits = factor(0:6)) # assign colors to values
      } else {
        gg <- ggplot(lls, aes(x = x, y = ll)) 
      }
      
      # for clipping the plot we need the extent
      min_ll <- min(lls$ll,ll_kappa)
      max_ll <- max(lls$ll,ll_kappa)
      
      gg  <- gg +  
        geom_line(linewidth = 1, group = 1) +
        geom_point(x = model_intercept, y = ll_hat, color = "red", size = 3) +
        geom_point(x = kappa_intercept, y = ll_kappa, color = "green", size = 5, shape = 1) +
        #    geom_vline(xintercept = mcmc_intercept, linetype = "solid", color = "darkblue") +
        ggtitle("optimal (red), sqrt(kappa) (green)") +
        theme_bw() +
        xlab(input$xaxis) +
        #  xlim(0,1) +
        theme( text=element_text(size=15),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
               legend.position = "none")  +
        ylab("log likelihood") +
        coord_cartesian(xlim = c(0,1), ylim = c(min_ll, max_ll), clip = 'off') 
      
      # if the bayesian data is available, shade the area under the curve
      
      return(gg)
    })
    
    
    # distribution for the mcmc
    output$mcmc_distro <- renderPlot({
      
      #input$action_tap # react to it
      input$action_tap

      req(mcmc_distro(), input$mcmc)
      if(!input$mcmc) return()
      
      my_stats <- mcmc() |> 
        filter(var == input$xaxis)
      
      
      my_title <- paste0("Mean: ", round(my_stats$avg,2),  
                         ", Mode 1: ", 
                         round(my_stats$mode1,2), 
                         ", Mode 2: ",round(my_stats$mode2,2) )
      
      pdf <-   mcmc_distro() |> 
        # renormalize log likelihood to be between 0 and 1
        mutate(likelihood = likelihood - min(likelihood),
               likelihood = likelihood / (max(likelihood) - min(likelihood))) |> 
        gather(var, value, -likelihood) |> 
        filter(var == input$xaxis)
      
      # average log likelihood plot data in bins of .05
      ll_pdf <- pdf |> 
        mutate(scaled = round(value*2, 1) * .5) |> 
        group_by(scaled) |> 
        summarise(likelihood = mean(likelihood),
                  value = mean(value)) |> 
        ungroup()
      
      pdf |> 
        ggplot(aes(x  = value, after_stat(scaled))) +
        geom_density(fill = "#8888FF") +
        geom_line(data = ll_pdf, aes(x = value, y = likelihood),  
                  color = "black", linetype = "dashed", group = 1) +
        theme_bw() +
        ggtitle(my_title) +
        theme(text=element_text(size=15),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        xlim(0,1)  +
        ylab("probability density") +
        coord_cartesian(clip = 'off')
      
      
    })
    
    
    
  })
}