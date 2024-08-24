# fleiss kappa
# Function to calculate Fleiss' kappa
fleiss_kappa <- function(counts) {
  
  matches <- function(N_c, N_r) {
    N_c_c  = N_r - N_c
    N_c*(N_c-1)/2 + N_c_c*(N_c_c-1)/2
  }
  
  # Number of cases (rows)
  n <- nrow(counts)
  # Number of raters (columns)
  k <- ncol(data)
  
  # Calculate the proportion of raters that assigned each category (0 or 1) to each case
  stats <- counts |> 
    mutate(p_ij = N_c / N_r,
           q_ij = 1 - p_ij,
           m = matches(N_c, N_r),
           p = m / matches(N_r,N_r)) |> 
    summarize(P_bar= mean(p),
              c = sum(N_c)/sum(N_r),
              P_e_bar= c^2 + (1-c)^2)
  
  # Calculate Fleiss' kappa
  kappa <- (stats$P_bar - stats$P_e_bar) / (1 - stats$P_e_bar)
  
  # see problems to kappa=0
  kappa <- if_else(kappa < 0 | is.na(kappa) | is.nan(kappa), 0, kappa)
  ll    <- -log_likelihood(c(stats$c, sqrt(kappa), stats$c), counts)
  degenerate <- is_degenerate(c(stats$c, sqrt(kappa), stats$c)) || kappa <= 0
  
  return(data.frame(t = stats$c, a = sqrt(kappa), p = stats$c, ll = ll, degenerate = degenerate))
}


# binomial mixture, returned as negative log likelihood
# negated because the optimizer minimizes
# N_r = number of raters
# N_c = number of raters who rated T_i = 1
# a = accuracy
# p = guess rate
fix_domain <- function(x){
  x[x < 0] <- 0
  x[x > 1] <- 1
  return(x)
}


negative_log_prob <- function(N_r, N_c, t, a, p){
  
  t <- fix_domain(t)
  a <- fix_domain(a)
  p <- fix_domain(p)
  
  prob_true <- dbinom(N_c, N_r, prob = a + (1-a)*p)
  prob_false <- dbinom(N_c, N_r, prob = (1-a)*p) # T = 0 case
  
  # despite the limits, the optimizer sometimes submits out-of-range parameters
  if(sum( is.nan(prob_true) + is.nan(prob_false)) > 0) return(1e6)
  
  # scale by t and (1-t) to get the mixture
  prob <- t*prob_true + (1-t)*prob_false
  prob[prob == 0] <- .00001
  
  return(-sum(log(prob), na.rm= TRUE)) # negated because the optimizer minimizes
}

# total log likelihood for a rating set
# x = parameters a and p
log_likelihood <- function(x, counts){
  t <- x[1]
  a <- x[2]
  p <- x[3]
  
  counts  %>%
    #rowwise() %>% # because log_prob isn't vectorized
    # note that the optimizer MINIMIZES
    summarize(ll = negative_log_prob(N_r, N_c, t, a, p)) %>%
    mutate(ll = if_else(is.infinite(ll), 1e6, ll)) %>%
    pull()
}

find_solution <- function(counts, init_params){
  
  params <- optim(par = init_params, 
                  fn = log_likelihood, 
                  counts = counts, 
                  lower = c(0,0,0), 
                  upper = c(1,1,1), 
                  method = "L-BFGS-B")$par
  
  return(c( params, -log_likelihood(params, counts)))
}

is_degenerate <- function(params, threshold = 1e-6){
  eps1 <- params[1]*params[2]*params[3]
  eps2 <- (1-params[1])*(1-params[2])*(1-params[3])
  
  if(eps1 > threshold & eps2 > threshold) {return(FALSE)}
  return(TRUE)
}

vector_line <- function(start, stop, steps){
  # create a line of points between two vectors
  return(tibble(t = seq(start[1], stop[1], length.out = steps),
                a = seq(start[2], stop[2], length.out = steps),
                p = seq(start[3], stop[3], length.out = steps)) |> 
           filter(row_number() > 1, row_number() < steps))
  
}

picket_fence <- function(counts){
  # create a grid of points around the edge of the parameter space
  # and return the negative log likelihood
  t <- c(.05, .1, .5, .9, .95)
  a <- c(.05, .1, .2)
  p <- c(.05, .1, .5, .9, .95)
  
  fence <- expand.grid(t = t, a = a, p = p)
  
  # apply the function log_likelihood to each row
  fence <- fence %>%
    rowwise() %>%
    mutate(ll = log_likelihood(c(t, a, p), counts)) |> 
    ungroup() |> 
    slice_min(ll, n = 1) |> 
    as.numeric()
  
  return(fence)
  
}

#' iterative optimization
iterative_optim <- function(counts, iterations = 1){
  # start with the middle case
  params <- find_solution(counts, c(.5,.5,.5))
  
  if(is_degenerate(params)) {
    fence_params <- picket_fence(counts)
    params <- find_solution(counts, fence_params)
  }
  
  return(data.frame( t = params[1], a = params[2], p = params[3], ll = params[4]))
  
}


# wrapper function to cache results
get_mcmc_results <- function(counts){
  hash_name <- rlang::hash(counts)
  # get filenames in the /cache folder
  cache_files <- list.files("cache/mcmc") |> 
    str_remove("\\.rds") 
  
  if(hash_name %in% cache_files){
    return(read_rds(paste0("cache/mcmc/",hash_name,".rds")))
  } else {
    results <- mcmc_stats(counts)
    write_rds(results, paste0("cache/mcmc/",hash_name,".rds"))
    return(results)
  }
  
}


# MCMC using stan
mcmc_stats <- function(counts){
  library(cmdstanr) 
  #  library(posterior)
  #  library(LaplacesDemon)
  
  compiled_model <- cmdstan_model("stan/t-a-p.stan")
  model_draws <- fit_tap_model(counts, compiled_model)
  tap_stats <- get_tap_stats(model_draws) |> 
    mutate(var = substr(var,1,1), # "accuracy" -> "a"
           Degenerate = !is.na(mode2)) 
  
  draws <- extract_vars(model_draws) |> 
    rename(a = accuracy)
  
  return(list(stats = tap_stats, draws = draws))
}

# MCMC for the t-a0,a1-p model
# wrapper function to cache results
get_mcmc_results2 <- function(counts, p0p1){
  hash_name <- rlang::hash(counts)
  # get filenames in the /cache folder
  if(isTRUE(p0p1)){
    folder <- "cache/mcmc3/"
  } else {
    folder <- "cache/mcmc2/"
  }
  
  cache_files <- list.files(folder) |> 
    str_remove("\\.rds") 
  
  if(hash_name %in% cache_files){
    return(read_rds(paste0(folder,hash_name,".rds")))
  } else {
    results <- mcmc_stats2(counts, p0p1)
    write_rds(results, paste0(folder,hash_name,".rds"))
    return(results)
  }
  
}


# MCMC using stan
mcmc_stats2 <- function(counts, p0p1){
  library(cmdstanr) 
  #  library(posterior)
  #  library(LaplacesDemon)
  
  if(isTRUE(p0p1)){
    stan_model <- "stan/t-a0a1-p0p1.stan"
  } else {
    stan_model <- "stan/t-a0a1-p.stan"
  }
  
  compiled_model <- cmdstan_model(stan_model)
  model_draws <- fit_tap_model(counts, compiled_model)
  tap_stats <- get_tap_stats(model_draws) |> 
    mutate(Degenerate = !is.na(mode2)) 
  
  draws <- extract_vars(model_draws) 
  
  return(list(stats = tap_stats, draws = draws))
}

#' Generate t-a-p statistics from a random effects data source
#' @param ratings a data frame that must include columns subject, rater, and C_ij (the binary classification)
#' @compiled_model a compiled stan model suitable for this data
#' @return a posterior stan results object
#' @details The model must be able to accomodate varying numbers of raters per subject
#' @export
fit_tap_model <- function(counts, compiled_model) {
  
  n_subjects <- n_distinct(counts$SubjectID__)
  
  N <- nrow(counts) # number of cases/subjects
  R <- counts$N_r   # ratings for each subject
  count <- counts$N_c # number of positive ratings
  
  fitted_model <- compiled_model$sample( 
    data = list(
      N = N, 
      R = R, 
      count = count), 
    seed = 123,
    chains = 3,
    parallel_chains = 3,
    refresh = 1000) 
  
  return(fitted_model$draws())
  
}

#' Auxillary function to extract the draws from a stan draws object
#' and put into a dataframe with columns for each variable
extract_vars <- function(model_draws){
  
  var_positions <- data.frame(position = 1:dim(model_draws)[3],
                              var = attr(model_draws, "dimnames")$variable) 
  
  # initialize a dataframe of the correct length by
  # storing the log probability column 
  
  draws_by_var <- tibble(likelihood = model_draws[,,1] |> as.vector() ) 
  
  for(i in 2:nrow(var_positions)){
    j <- var_positions$position[i]
    var_name <- var_positions$var[i]
    
    draws_by_var <- draws_by_var %>%
      mutate(!!var_name := (model_draws[,,j] |> as.vector()))
  }
  
  return(draws_by_var)
}

#' Get t-a-p coefficients
#' @param compiled_model The stan model for t-a0,a1-p_i
#' @param model_data A dataframe with columns: ID n_ratings, and sum_ratings. 
#' ID is the subject identifier, n_ratings is the number of binary ratings,
#' and sum_ratings is the sum of the binary ratings.
#' @return a list with the t-a-p coefficients. If there are two modes for a
#' these are returned as a0 and a1
get_tap_stats <- function(model_draws){
  
  model_means <- data.frame(position = 1:dim(model_draws)[3],
                            var = attr(model_draws, "dimnames")$variable,
                            avg = NA, # mean 
                            p05 = NA, # 5th percentile
                            p25 = NA, # 25th percentile
                            median = NA, # median
                            p75 = NA, # 75th percentile
                            p95 = NA, # 95th percentile
                            mode1 = NA,
                            mode2 = NA) |> 
    filter(str_detect(var,"^t|^a|^p")) 
  
  # get parameter averages
  for(i in 1:nrow(model_means)){
    j <- model_means$position[i]
    
    model_means$avg[i] <- model_draws[,,j] %>% mean()
    model_means$median[i] <- model_draws[,,j] %>% median()
    model_means$p05[i] <- quantile(model_draws[,,j], 0.05)
    model_means$p25[i] <- quantile(model_draws[,,j], 0.25)
    model_means$p75[i] <- quantile(model_draws[,,j], 0.75)
    model_means$p95[i] <- quantile(model_draws[,,j], 0.95)
    
    # get modes
    model_means$mode1[i] <- Modes(model_draws[,,j])$modes[1]
    
    if(length(Modes(model_draws[,,j])$modes) > 1){
      model_means$mode2[i] <- Modes(model_draws[,,j])$modes[2]
    }
  }
  
  return(model_means)
}

# binomial mixture, returned as negative log likelihood
# negated because the optimizer minimizes
# N_r = number of raters
# N_c = number of raters who rated T_i = 1
# a = accuracy
# p = guess rate
fix_domain <- function(x){
  x[x < 0] <- 0
  x[x > 1] <- 1
  return(x)
}


negative_log_prob <- function(N_r, N_c, t, a, p){
  
  t <- fix_domain(t)
  a <- fix_domain(a)
  p <- fix_domain(p)
  
  prob_true <- dbinom(N_c, N_r, prob = a + (1-a)*p)
  prob_false <- dbinom(N_c, N_r, prob = (1-a)*p) # T = 0 case
  
  # despite the limits, the optimizer sometimes submits out-of-range parameters
  if(sum( is.nan(prob_true) + is.nan(prob_false)) > 0) return(1e6)
  
  # scale by t and (1-t) to get the mixture
  prob <- t*prob_true + (1-t)*prob_false
  prob[prob == 0] <- .00001
  
  return(-sum(log(prob), na.rm= TRUE)) # negated because the optimizer minimizes
}

# total log likelihood for a rating set
# x = parameters a and p
log_likelihood <- function(x, counts){
  t <- x[1]
  a <- x[2]
  p <- x[3]
  
  counts  %>%
    #rowwise() %>% # because log_prob isn't vectorized
    # note that the optimizer MINIMIZES
    summarize(ll = negative_log_prob(N_r, N_c, t, a, p)) %>%
    mutate(ll = if_else(is.infinite(ll), 1e6, ll)) %>%
    pull()
}

find_solution <- function(counts, init_params){
  
  params <- optim(par = init_params, 
                  fn = log_likelihood, 
                  counts = counts, 
                  lower = c(0,0,0), 
                  upper = c(1,1,1), 
                  method = "L-BFGS-B")$par
  
  return(c( params, log_likelihood(params, counts)))
}

is_degenerate <- function(params, threshold = 1e-6){
  # check for blanks
  if(sum(is.na(params)) > 0) return(TRUE)
  
  eps1 <- params[1]*params[2]*params[3]
  eps2 <- (1-params[1])*(1-params[2])*(1-params[3])
  
  if(eps1 > threshold & eps2 > threshold) {return(FALSE)}
  return(TRUE)
}

vector_line <- function(start, stop, steps){
  # create a line of points between two vectors
  return(tibble(t = seq(start[1], stop[1], length.out = steps),
                a = seq(start[2], stop[2], length.out = steps),
                p = seq(start[3], stop[3], length.out = steps)) |> 
           filter(row_number() > 1, row_number() < steps))
  
}

picket_fence <- function(counts){
  # create a grid of points around the edge of the parameter space
  # and return the negative log likelihood
  t <- c(.05, .1, .5, .9, .95)
  a <- c(.05, .1, .2)
  p <- c(.05, .1, .5, .9, .95)
  
  fence <- expand.grid(t = t, a = a, p = p)
  
  # apply the function log_likelihood to each row
  fence <- fence %>%
    rowwise() %>%
    mutate(ll = log_likelihood(c(t, a, p), counts)) |> 
    ungroup() |> 
    slice_min(ll, n = 1) |> 
    as.numeric()
  
  return(fence)
  
}

#' iterative optimization
iterative_optim <- function(counts, iterations = 1){
  # start with the middle case
  params <- find_solution(counts, c(.5,.5,.5))
  
  if(is_degenerate(params)) {
    fence_params <- picket_fence(counts)[1:3]
    params <- find_solution(counts, fence_params)
  }
  
  return(data.frame( t = params[1], 
                     a = params[2], 
                     p = params[3], 
                     ll = -params[4],
                     degenerate = is_degenerate(params)))
  
}

#' faceted density plots for mcmc draws
get_density <- function(x){
  d <- density(x, from = 0, to = 1, bw = .02)
  q <- quantile(x, c(0.05, 0.25, 0.45, .55, 0.75, 0.95)) |> as.list()
  
  out <- data.frame(x = d$x, y = d$y) |> 
    mutate(y = y/max(y), # scale to 0,1
           region = case_when(x < q$`5%` ~ 0,
                              x < q$`25%` ~ 1,
                              x < q$`45%` ~ 2,
                              x < q$`55%` ~ 3,
                              x < q$`75%` ~ 4,
                              x < q$`95%`~ 5, 
                              TRUE ~ 6),
           region = as.factor(region))
  
  return(out)
}

plot_draw_densities <- function(draws){
  #  draws <- read_rds("shiny/cache/mcmc3/727e989915209e1d9ad556facca13c41.rds")$draws
  
  # prepare for plotting
  pdf <- draws |> 
    gather(var, value, -likelihood) 
  
  # change title based on the vars present
  if(any(str_detect(pdf$var,"p0"))) {
    my_title <- "t-a0,a1-p0,p1"
  } else {
    my_title <- "t-a0,a1-p"
  }
  
  
  param_density <- pdf |> 
    group_by(var) %>%
    do(get_density(.$value))
  
  # create likelihood densities
  ll_density <- pdf |>
    group_by(var) |> 
    mutate(scaled = round(value*2, 1) * .5) |> 
    group_by(var, scaled) |> 
    summarise(likelihood = mean(likelihood),
              value = mean(value)) |> 
    group_by(var) |> 
    mutate(likelihood = (likelihood - min(likelihood))/(max(likelihood) - min(likelihood)) ) |> 
    ungroup() |> 
    select(-scaled)
  
  ll_means <- pdf |> 
    group_by(var) |> 
    summarise(value = mean(value))
  
  param_density |> 
    ggplot(aes(x  = x, y = y)) +
    geom_line()  +
    geom_ribbon( aes( fill = region, ymax = y, group = region), 
                 ymin = 0, 
                 position = position_identity()) +
    scale_fill_manual(values = c("#EEEEFF33", 
                                 "#BBBBFF55", 
                                 "#4444FF88", 
                                 "#2222AA",
                                 "#4444FF88", 
                                 "#BBBBFF55", 
                                 "#EEEEFF33"),
                      drop = TRUE,  # omit unused factors
                      limits = factor(0:6)) +
    geom_line(data = ll_density, aes(x = value, y = likelihood), 
              linetype = "dotted", color = "#666666") +
    geom_vline(data = ll_means, aes(xintercept = value), color = "darkorange") +
    geom_label(data = ll_means, 
               aes(x = value, label = round(value,2)), 
               y = .5, color = "orange",
               label.padding = unit(0.15, "lines")) +
    ggtitle(my_title) +
    theme_bw() +
    theme(text=element_text(size=15),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none") +
    xlim(0,1)  +
    ylab("probability density") +
    xlab("") +
    facet_grid(var ~ ., scales = "free_y")
}

#' tap optim over ordinal scale
#' @param ratings the standarized long format ratings df 
#'                with subject id and a rating column 
get_ordinal_tap <- function(ratings){
  rating_values <- sort(unique(ratings$rating))
  n_values <- length(rating_values)
  
  output <- data.frame()
  for(i in 1:(n_values - 1)){
    lower_rating <- rating_values[i]
    upper_rating <- rating_values[i+1]
    
    counts <- ratings |> 
      group_by(SubjectID__) |>
      summarize(N_r = n(),
                N_c = sum(rating <= lower_rating)) |> 
      filter(N_r >= 2)
    
    tap_params <- iterative_optim(counts) |> 
      mutate(CutPoint = str_c(lower_rating,"|", upper_rating),
             type = "t-a-p")
    
    fleiss_params <- fleiss_kappa(counts) |> 
      mutate(CutPoint = str_c(lower_rating,"|", upper_rating),
             type = "Fleiss")
    
    output <- output |> 
      rbind(tap_params, fleiss_params)
    
  }
  
  return(output)
}