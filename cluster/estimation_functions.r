library(sl3)
library(MBESS)

# define a function to estimate the regular cohen's d
estimate_cohens_d <- function(data) {
    # estimate the cohen's d
    e_y_1 <- mean(data$y[data$a == 1])
    e_y_0 <- mean(data$y[data$a == 0])
    v_y_1 <- var(data$y[data$a == 1])
    v_y_0 <- var(data$y[data$a == 0])
    pooled_v <-
    ((sum(data$a == 1) - 1) * v_y_1 + (sum(data$a == 0) - 1) * v_y_0) /
    (nrow(data) - 2)
    cohens_d <- (e_y_1 - e_y_0) / sqrt(pooled_v)

    # calculate the confidence interval
    ci_result <- ci.smd(
        smd = cohens_d, n.1 = sum(data$a == 1), n.2 = sum(data$a == 0)
        )
    cohens_d_lb <- ci_result$Lower.Conf.Limit.smd
    cohens_d_ub <- ci_result$Upper.Conf.Limit.smd

    return(data.frame(
        cohens_d = cohens_d,
        cohens_d_lb = cohens_d_lb,
        cohens_d_ub = cohens_d_ub
    ))
}

# define a function to estimate the corresponding causal ES
# use superlearners to estimate the nuisance parameters
estimate_causal_es <- function(data) {
   # define individual learners
   lrnr_rf <- Lrnr_ranger$new(mtry = 2)
   lrnr_glm <- Lrnr_glm$new()
   lrnr_gam <- Lrnr_gam$new()
   lrnr_xgb <- Lrnr_xgboost$new()

   # stack the learners
  stack <- Stack$new(lrnr_rf, lrnr_glm, lrnr_gam, lrnr_xgb)

   # make a sl for the propensity score
   sl_ps <- Lrnr_sl$new(
    learners = stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_loglik_binomial)
   )

   # define the task for propensity score estimation
   task_ps <- make_sl3_Task(
    data = data,
    outcome = "a",
    covariates = c("w1", "w2")
   )

   # make a sl for regression taks
   sl_reg <- Lrnr_sl$new(
    learners = stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_squared_error)
   )

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

   # train the superlearners
   fit_ps <- sl_ps$train(task_ps)
   fit_reg <- sl_reg$train(task_outcome)
   fit_regsq <- sl_reg$train(task_outcomesq)

   # get the propensity score for each observation
   # this is for the influence curve
   ps <- fit_ps$predict(task_ps)

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
    
    # Restore original data for influence curve calculations
    data$a <- a_original

    # cauculate the plug-in estimate
    v_y_1 <- mean(q_1_sq) - mean(q_1)^2
    v_y_0 <- mean(q_0_sq) - mean(q_0)^2
    g_1 <- sum(a_original) / nrow(data)
    g_0 <- 1 - g_1
    es_plugin <- (mean(q_1) - mean(q_0)) / sqrt(g_0 * v_y_0 + g_1 * v_y_1)
    
    # estimate the influence curve
    d_1_0 <- (a_original == 0) / (1 - ps) * (data$y - q_0) + q_0 - mean(q_0)
    d_1_1 <- (a_original == 1) / ps * (data$y - q_1) + q_1 - mean(q_1)
    d_2_0 <- (a_original == 0) / (1 - ps) * (data$y_sq - q_0_sq) + q_0_sq - mean(q_0_sq)
    d_2_1 <- (a_original == 1) / ps * (data$y_sq - q_1_sq) + q_1_sq - mean(q_1_sq)
    d_g_0 <- (a_original == 0) - g_0

    # calculate the correction term
    z <- g_0 * v_y_0 + g_1 * v_y_1
    ic_1_1 <- (g_0 * v_y_0 + g_1 * (mean(q_1_sq) - mean(q_1) * mean(q_0))) / z ^ (1.5) * d_1_1
    ic_1_0 <- ((-g_0) * (mean(q_0_sq) - mean(q_0) * mean(q_1)) - g_1 * v_y_1) / z ^ (1.5) * d_1_0
    ic_2_1 <- g_1 * (mean(q_1) - mean(q_0)) / z ^ (1.5) * d_2_1 * (-1/2)
    ic_2_0 <- g_0 * (mean(q_1) - mean(q_0)) / z ^ (1.5) * d_2_0 * (-1/2)
    ic_g_0 <- (mean(q_1) - mean(q_0)) * (v_y_0 - v_y_1) / z ^ (1.5) * d_g_0 * (-1/2)
    ic <- ic_1_0 + ic_1_1 + ic_2_0 + ic_2_1 + ic_g_0

    # calculate the one-step estimate
    es_one_step <- es_plugin + mean(ic)

    # calculate the confidence interval
    ic_var <- var(ic)
    ic_se <- sqrt(ic_var / nrow(data))
    es_one_step_lb <- es_one_step - 1.96 * ic_se
    es_one_step_ub <- es_one_step + 1.96 * ic_se

    return(data.frame(
        es_plugin = es_plugin,
        es_one_step = es_one_step,
        es_one_step_lb = es_one_step_lb,
        es_one_step_ub = es_one_step_ub
    ))
}

estimate_resi <- function(data){
  
}