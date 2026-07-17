library(sl3)
library(MBESS)
library(RESI)

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

# Estimate the equal-variance-weighted causal effect size and causal ATE RESI.
# Nuisance functions are estimated with outer two-fold cross-fitting.
estimate_ces_and_resi <- function(data, ps_bound = 0.01, conf_level = 0.95) {
  n <- nrow(data)
  fold_id <- integer(n)
  for (treatment in 0:1) {
    indices <- which(data$a == treatment)
    fold_id[indices] <- sample(rep(1:2, length.out = length(indices)))
  }

  ps <- q_1 <- q_0 <- q_1_sq <- q_0_sq <- rep(NA_real_, n)

  make_stack <- function() {
    Stack$new(
      Lrnr_ranger$new(mtry = 2, num.threads = 1),
      Lrnr_glm$new(),
      Lrnr_gam$new(),
      Lrnr_xgboost$new(nthread = 1)
    )
  }

  make_prediction_task <- function(prediction_data, outcome) {
    make_sl3_Task(
      data = prediction_data,
      outcome = outcome,
      covariates = c("w1", "w2", "a")
    )
  }

  for (fold in 1:2) {
    validation_indices <- which(fold_id == fold)
    training_data <- data[fold_id != fold, , drop = FALSE]
    validation_data <- data[validation_indices, , drop = FALSE]

    ps_learner <- Lrnr_sl$new(
      learners = make_stack(),
      metalearner = Lrnr_nnls$new(eval_function = loss_loglik_binomial)
    )
    outcome_learner <- Lrnr_sl$new(
      learners = make_stack(),
      metalearner = Lrnr_nnls$new(eval_function = loss_squared_error)
    )
    outcome_sq_learner <- Lrnr_sl$new(
      learners = make_stack(),
      metalearner = Lrnr_nnls$new(eval_function = loss_squared_error)
    )

    ps_training_task <- make_sl3_Task(
      data = training_data,
      outcome = "a",
      covariates = c("w1", "w2")
    )
    outcome_training_task <- make_prediction_task(training_data, "y")
    outcome_sq_training_task <- make_prediction_task(training_data, "y_sq")

    ps_fit <- ps_learner$train(ps_training_task)
    outcome_fit <- outcome_learner$train(outcome_training_task)
    outcome_sq_fit <- outcome_sq_learner$train(outcome_sq_training_task)

    ps_validation_task <- make_sl3_Task(
      data = validation_data,
      outcome = "a",
      covariates = c("w1", "w2")
    )
    ps[validation_indices] <- ps_fit$predict(ps_validation_task)

    validation_data_1 <- validation_data
    validation_data_1$a <- 1
    validation_data_0 <- validation_data
    validation_data_0$a <- 0

    q_1[validation_indices] <- outcome_fit$predict(
      make_prediction_task(validation_data_1, "y")
    )
    q_0[validation_indices] <- outcome_fit$predict(
      make_prediction_task(validation_data_0, "y")
    )
    q_1_sq[validation_indices] <- outcome_sq_fit$predict(
      make_prediction_task(validation_data_1, "y_sq")
    )
    q_0_sq[validation_indices] <- outcome_sq_fit$predict(
      make_prediction_task(validation_data_0, "y_sq")
    )
  }


  ps <- pmin(pmax(ps, ps_bound), 1 - ps_bound)

  y_1_os <- (data$a == 1) / ps * (data$y - q_1) + q_1
  y_0_os <- (data$a == 0) / (1 - ps) * (data$y - q_0) + q_0
  y_1_sq_os <- (data$a == 1) / ps * (data$y_sq - q_1_sq) + q_1_sq
  y_0_sq_os <- (data$a == 0) / (1 - ps) * (data$y_sq - q_0_sq) + q_0_sq

  mu_1 <- mean(y_1_os)
  mu_0 <- mean(y_0_os)
  nu_1 <- mean(y_1_sq_os)
  nu_0 <- mean(y_0_sq_os)
  v_y_1 <- nu_1 - mu_1 ^ 2
  v_y_0 <- nu_0 - mu_0 ^ 2
  g_1 <- g_0 <- 0.5
  z <- g_0 * v_y_0 + g_1 * v_y_1
  if (!is.finite(z) || z <= 0) {
    stop("estimated pooled potential-outcome variance is not positive")
  }

  ate <- mu_1 - mu_0
  est <- ate / sqrt(z)

  d_mu_1 <- y_1_os - mu_1
  d_mu_0 <- y_0_os - mu_0
  d_nu_1 <- y_1_sq_os - nu_1
  d_nu_0 <- y_0_sq_os - nu_0
  ic <- (d_mu_1 - d_mu_0) / sqrt(z) -
    ate / (2 * z ^ 1.5) * (
      g_1 * (d_nu_1 - 2 * mu_1 * d_mu_1) +
        g_0 * (d_nu_0 - 2 * mu_0 * d_mu_0)
    )

  alpha <- 1 - conf_level
  critical_value <- qnorm(1 - alpha / 2)
  est_se <- sqrt(var(ic) / n)
  est_lb <- est - critical_value * est_se
  est_ub <- est + critical_value * est_se

  ate_ic <- y_1_os - y_0_os - ate
  resi_result <- resi_estimation(
    ate = ate,
    ate_ic = ate_ic,
    n = n,
    conf_level = conf_level
  )

  data.frame(
    est = est,
    est_se = est_se,
    est_lb = est_lb,
    est_ub = est_ub,
    resi = resi_result$resi,
    resi_lb = resi_result$resi_lb,
    resi_ub = resi_result$resi_ub
  )
}

# Bias-corrected scalar RESI and a noncentral chi-square confidence interval.
resi_estimation <- function(ate, ate_ic, n, conf_level = 0.95) {
  ate_ic_variance <- var(ate_ic)
  if (!is.finite(ate_ic_variance) || ate_ic_variance <= 0) {
    stop("the estimated ATE influence-curve variance is not positive")
  }

  wald_statistic <- n * ate ^ 2 / ate_ic_variance
  resi <- sqrt(max(0, (wald_statistic - 1) / n))
  ncp_limits <- suppressWarnings(
    conf.limits.nc.chisq(
      Chi.Square = wald_statistic,
      conf.level = conf_level,
      df = 1
    )
  )

  data.frame(
    resi = resi,
    resi_lb = sqrt(ncp_limits$Lower.Limit / n),
    resi_ub = sqrt(ncp_limits$Upper.Limit / n)
  )
}