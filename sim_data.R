data_sim_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel( 

    sliderInput(ns("N_s"),
                "Number of subjects",
                min = 2,
                max = 1000,
                value = 100),
    sliderInput(ns("N_r"),
                "Number of raters",
                min = 2,
                max = 100,
                value = 10),
    # t parameter
    sliderInput(ns("t"),
                "True rate of Class 1 (t)",
                min = 0,
                max = 1,
                value = .5),
    
    # a parameter
    checkboxInput(ns("a_type"), "Use a0,a1", value = FALSE),
    uiOutput(ns( "slider_a")),
    
    
    # p parameter
    checkboxInput(ns("p_proficient"), "Unbiased raters (p = t)", value = FALSE),
    uiOutput(ns( "nested_p_ui")),
    #checkboxInput(ns("p_type"), "Use p0,p1", value = FALSE),
    
    
    # should we produce the exact distribution or sample from it?
    radioButtons(ns("sample"), "Data type", c("Exact", "Random"), inline = TRUE),
    actionButton(ns( "action"), "Generate data")),
    mainPanel(
        tableOutput(ns("match_stats")),
        plotOutput(ns("histogram")),
        textOutput(ns("description"))
    )
  )
}

data_sim_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    sim_ratings <- reactiveVal(NULL)
    
    # UI portion -- dynamic slider inputs
    output$slider_a <- renderUI({
      if (input$a_type) {
        tagList(
          sliderInput(ns("a0"), "Accuracy for Class 0 (a0)", min = 0, max = 1, value = .5),
          sliderInput(ns( "a1"), "Accuracy for Class 1 (a1)", min = 0, max = 1, value = .5)
        )
      } else {
        sliderInput(ns( "a"), "Accuracy (a)", min = 0, max = 1, value = .5)
      }
    })
    
    
    # should we set p = c? If so we don't need the p sliders
    output$nested_p_ui <- renderUI({
      if (!input$p_proficient) {
        tagList(
          checkboxInput(ns( "p_type"), "Use p0,p1", value = FALSE),
          uiOutput(ns( "p_sliders"))
        )
      }
    })
    
    # p slider or p0,p1 sliders, depending on the checkbox
    output$p_sliders <- renderUI({
      
      if (input$p_type) {
        tagList(
          sliderInput(ns("p0"), "Chance of Class 0 (p0)", min = 0, max = 1, value = .5),
          sliderInput(ns( "p1"), "Chance of Class 1 (p1)", min = 0, max = 1, value = .5)
        )
      } else {
        sliderInput(ns( "p"), "Chance of Class 1 (p)", min = 0, max = 1, value = .5)
      }
      
    })
    
    # print out the match statistics for the slider settings
    output$match_stats <- renderTable({
      params <- get_params(input)
      
      # convenience function for formatting output
      fmt <- function(x, y= NULL, ratio = FALSE) {
        if(is.null(y)) {
          sprintf("%.2f", x)
        } else if(!ratio) {
          sprintf("%.2f / %.2f", x, y)
        } else{
          sprintf("%.2f / %.2f = %.2f", x, y, x/y)
        }
        
      }
      
      with(params, 
           tibble(`Class 1/0 rating proportion` = fmt( t*(a1 + (1-a1)*p1) + (1-t)*((1-a0)*p0),
                                                                        1 - (t*(a1 + (1-a1)*p1) + (1-t)*((1-a0)*p0)) ),
                  `Accurate match rate of Class 1/0` =  fmt(t*a1*a1, (1-t)*a0*a0),
                  `Total accurate match rate` = fmt(t*a1*a1 + (1-t)*a0*a0),
                  `Random match rate of Class 1/0` = fmt(t*((1-a1)*p1)^2 + (1-t)*((1-a0)*p0)^2,
                                                         t*((1-a1)*(1-p1))^2 + (1-t)*((1-a0)*(1-p0))^2),
                  `Mixed match rate of Class 1/0` = fmt(2*t*a1*(1-a1)*p1, # subject is class 1
                                                        2*(1-t)*a0*(1-a0)*(1-p0)), # subject is class 0
                  `Accurate / total match rate` = fmt(t*a1*a1 + (1-t)*a0*a0,
                                                      t*a1*a1 + (1-t)*a0*a0 +
                                                      t*((1-a1)*p1)^2 + (1-t)*((1-a0)*p0)^2 + 
                                                      t*((1-a1)*(1-p1))^2 + (1-t)*((1-a0)*(1-p0))^2 +
                                                      2*t*a1*(1-a1)*p1 +
                                                      2*(1-t)*a0*(1-a0)*(1-p0), TRUE)
                   ) |> 
             gather(`Asymptotic statistic`, `Value`) 
      )
    })
 
    
    
    # Take a reactive dependency on input$action, but
    # not on any of the stuff inside the function
    observeEvent(input$action,{
      
      # get count params
      N_r = isolate(input$N_r)
      N_s = isolate(input$N_s)
      
      # create short names for the parameters
      params <- get_params(isolate(input))
      
      with(params,
        if(input$sample == "Exact") {
          tdf <<- generate_exact_ratings(N_s, N_r, t, a0, a1, p0, p1)
        } else {
          tdf <<- generate_sample_ratings(N_s, N_r, t, a0, a1, p0, p1)
        }
      )
      
      # make all columns character vars
      tdf <- tdf |> 
            mutate(across(everything(), as.character))
      
      # set the reactive value to df, which should trigger events
      sim_ratings(tdf)
    }) 
      
    output$histogram <- renderPlot({
        req(sim_ratings())
      
      # get count params
      N_r = isolate(input$N_r)
      N_s = isolate(input$N_s)
      
      # create short names for the parameters
      params <- get_params(isolate(input))
      
      # calculate the stats we need
      with(params,{
        c0_probs <<- dbinom(0:N_r, N_r, prob = (1-a0)*p0)
        c1_probs <<- dbinom(0:N_r, N_r, prob = a1 + (1-a1)*p1)
        
        c_0_mean <<- N_r*(1-a0)*p0
        c_1_mean <<- N_r*(a1 + (1-a1)*p1)

        tap_prob <<- c0_probs*(1-t) + c1_probs*t
      }
      )

      rating_sum_counts <- data.frame(count_1s = 0:N_r, 
                                      y = round_preserve_sum(tap_prob*N_s))
        
      # Calculate counts of 1 ratings per subject
      rating_counts <- sim_ratings() %>%
          group_by(SubjectID__) %>%
          summarize(count_1s = sum(rating == 1))
        
      # Create the histogram
      ggplot(rating_counts, aes(x = count_1s)) +
          geom_histogram(binwidth = 1, color = "white", fill = "steelblue") +
          geom_line(data = rating_sum_counts, aes(x = count_1s, y = y), 
                    color = "darkorange", linewidth = 1.5, linetype = "dashed") +
          geom_point(x = c_0_mean, y = 0, shape = 18, color = "black", size = 4) +
          geom_point(x = c_1_mean, y = 0, shape = 18, color = "gray", size = 4) +
          labs(title = "Histogram of Class 1 Ratings per Subject", 
               x = "Count of Class 1 Ratings", y = "Frequency") +
          theme_bw()
      })
    
    output$description <- renderText({
      req(sim_ratings())
      "The histogram shows the per-subject counts of Class 1 ratings from the
          generated data, with the dashed line representing the expected distribution. 
          The gray and black markers show the expected means for Class 1 
          and Class 0, repectively. If the black marker is to the right
          of the gray marker, it means the model is subject to label-switching,
          where we can't tell which of the two categories is which from
          the data."
    })
    
    return(sim_ratings)
  })
}


#' simplify parameters to t, a0,a1,p0,p1
#' @param input
#' @return a list of simplified parameters
get_params <- function(input){
  # make sure checkboxes are initialized
  if(is.null(input$p_proficient)) return()
     
  t = input$t
  
  if(input$a_type) {
    a0 <- input$a0
    a1 <- input$a1
  } else {
    a0 <- input$a
    a1 <- input$a
  }
  
  if(input$p_proficient) {
    p0 <- p1 <- t
    
  } else {
  
    if(isTruthy(input$p_type)) { # it might be null
      p0 <- input$p0
      p1 <- input$p1
    } else {
      p0 <- input$p
      p1 <- input$p
    }
  }
  
  list(N_r = input$N_r, N_s = input$N_s, t = t, a0 = a0, a1 = a1, p0 = p0, p1 = p1)
}

#' get ratings by sampling a known distribution
#' @param n_subjects number of subjects
#' @param n_raters number of raters per subject on average
#' @param t truth average
#' @param a accuracy average
#' @param p guess rate average
#' @return a data frame with ratings and parameters, including the input t-a-p averages, the 
#' random effect parameters p_i and a_j and the discrete random variables T_i, P_ij, and A_ij. 
#' @details This function generates ratings for a set of subjects based on a t-a-p
#' model. It uses t to determine the 0 or 1 truth values for each subject, 
#' creates random effects p_i and a_j for each subject and rater, respectively, using a beta distribution.
#' @export
generate_sample_ratings <- function(n_subjects = 100, n_raters = 5, 
                                    t = .5, a0 = .7, a1 = .7, 
                                    p0 = .5, p1 = .5) {
  
  subject_params <- tibble(SubjectID__ = 1:n_subjects,   
                           T_i = sample(0:1,
                                        n_subjects, 
                                        replace = TRUE, 
                                        prob = c(1 - t, t))) 
  

  rater_params <- tibble(rater = 1:n_raters)
 
  
  param_grid <- subject_params %>%
    cross_join(rater_params)
  
  n_ratings <- nrow(param_grid)
  
  # generate ratings based on t_i-a_j-p_i model
  param_grid <- param_grid %>% 
    mutate(a_j = if_else(T_i == 1, a1, a0),
           p_i = if_else(T_i == 1, p1, p0),
           A_ij = as.integer(a_j > runif(n_ratings)),
           P_ij = as.integer(p_i > runif(n_ratings)),
           C_ij = if_else(A_ij == 1,T_i,P_ij)) |> 
    # reorder output
    select(SubjectID__, rating = C_ij) # random variables
  
  return(param_grid)
}

# cf https://www.r-bloggers.com/2016/07/round-values-while-preserve-their-rounded-sum-in-r/
round_preserve_sum <- function(x, digits = 0) {
  up = 10 ^ digits
  x = x * up
  y = floor(x)
  indices = tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] = y[indices] + 1
  y / up
}


#' get ratings that have a given distribution
#' @param n_subjects number of subjects
#' @param n_raters number of raters per subject on average
#' @param t truth average
#' @param a accuracy average
#' @param p guess rate average
#' @return a data frame with ratings and parameters, including the input t-a-p averages, the 
#' random effect parameters p_i and a_j and the discrete random variables T_i, P_ij, and A_ij. 
#' @details This function generates ratings for a set of subjects based on a t-a-p
#' model. It uses t to determine the 0 or 1 truth values for each subject, 
#' creates random effects p_i and a_j for each subject and rater, respectively, using a beta distribution.
#' @export
generate_exact_ratings <- function(N_s, N_r, t, a0, a1, p0, p1) {
  # how many values are in the distro?
  
  # cumulative probabilities as cut points, with leading zero
  c0_probs <- dbinom(0:N_r, N_r, prob = (1-a0)*p0)
  c1_probs <- dbinom(0:N_r, N_r, prob = a1 + (1-a1)*p1) 
  
  tap_prob <- c0_probs*(1-t) + c1_probs*t
  
  rating_sum_counts <- round_preserve_sum(tap_prob*N_s)
  
  # now that we have counts, produce the ratings
  
  param_grid <- expand.grid(rater = 1:N_r,SubjectID__ = 1:N_s ) 
  
  ratings <- c()
  
  for(i in 0:(length(rating_sum_counts)-1)) {
    
    # make a single rating of this type              
    r <- c(rep(1, i), rep(0, N_r - i))
    
    # repeat that rating for the number of subjects
    r <- rep(r, rating_sum_counts[i + 1])
    
    #add it to the list
    ratings <- c(ratings, r)
  }
  
  param_grid <- param_grid |> 
    mutate(rating = as.integer(ratings)) |> 
    select(-rater)
  
  return(param_grid)
}
  
  