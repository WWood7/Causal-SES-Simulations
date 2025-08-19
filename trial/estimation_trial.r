library(sl3)


# estimate the plug-in estimates using rf learners
# investigate the effect of mtry and depth
estimate_causal_es <- function(data) {
   # define the regression tasks
   # the outcome regression
   task_outcome <- make_sl3_Task(
    data = data,
    outcome = "y",
    covariates = c("w1", "w2", "a")
   )

   # the squared outcome regression
   task_outcomesq <- make_sl3_Task(
    data = data,
    outcome = "y_sq",
    covariates = c("w1", "w2", "a")
   )

  # now define a grid of hyperparameters for the random forest learners
  mtry_grid <- c(1, 2)
  depth_grid <- c(2, 5, 10, 20, 50)
  grid <- expand.grid(mtry = mtry_grid, depth = depth_grid)

  # initialize a list to store the results
  results <- list()

  # now loop through the grid and train the random forest learners
  for (i in 1:nrow(grid)) {
    mtry <- grid$mtry[i]
    depth <- grid$depth[i]
    lrnr_rf <- Lrnr_ranger$new(mtry = mtry, max_depth = depth)
    fit_reg <- lrnr_rf$train(task_outcome)
    fit_regsq <- lrnr_rf$train(task_outcomesq)

    # calculate the plug-in estimates
    # conserve the original A
    a_original <- data$a

    # Create counterfactual datasets
    data_1 <- data
    data_1$a <- 1
    data_0 <- data
    data_0$a <- 0
    
    # Create new tasks for counterfactual predictions
    task_cf_1 <- make_sl3_Task(
      data = data_1,
      outcome = "y",
      covariates = c("w1", "w2", "a")
    )
    
    task_cf_1_sq <- make_sl3_Task(
      data = data_1,
      outcome = "y_sq",
      covariates = c("w1", "w2", "a")
    )
    
    task_cf_0 <- make_sl3_Task(
      data = data_0,
      outcome = "y",
      covariates = c("w1", "w2", "a")
    )
    
    task_cf_0_sq <- make_sl3_Task(
      data = data_0,
      outcome = "y_sq",
      covariates = c("w1", "w2", "a")
    )

    # outcome regression with a = 1
    q_1 <- fit_reg$predict(task_cf_1)
    q_1_sq <- fit_regsq$predict(task_cf_1_sq)

    # outcome regression with a = 0
    q_0 <- fit_reg$predict(task_cf_0)
    q_0_sq <- fit_regsq$predict(task_cf_0_sq)
    
    # calculate the plug-in estimate
    v_y_1 <- mean(q_1_sq) - mean(q_1)^2
    v_y_0 <- mean(q_0_sq) - mean(q_0)^2
    g_1 <- sum(a_original) / nrow(data)
    g_0 <- 1 - g_1
    es_plugin <- (mean(q_1) - mean(q_0)) / sqrt(g_0 * v_y_0 + g_1 * v_y_1)

    # store the results
    results[[i]] <- data.frame(
      mtry = mtry,
      depth = depth,
      bias = es_plugin - 0.876
    )
  }
  # merge the results
  results <- do.call(rbind, results)
  return(results)
}