# this script is used to check the performance of the ranger learner
library(ranger)
library(sl3)


source("data_generation.r")

data <- generate_data(n = 2000, effect_size_type = 3)

# define the learners
lrnr_rf <- Lrnr_ranger$new(num.trees = 1000)
lrnr_glm <- Lrnr_glm$new()
lrnr_gam <- Lrnr_gam$new()
lrnr_xgb <- Lrnr_xgboost$new()
lrnr_hal <- Lrnr_hal9001$new()

# stack the learners
stack <- Stack$new(lrnr_rf, lrnr_glm, lrnr_gam, lrnr_xgb, lrnr_hal)

# make a sl for propensity score estimation
sl_ps <- Lrnr_sl$new(
    learners = stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_loglik_binomial)
)

# make a sl for regression
sl_reg <- Lrnr_sl$new(
    learners = stack,
    metalearner = Lrnr_nnls$new(eval_function = loss_squared_error)
)

# define the task for regression
task_reg <- make_sl3_Task(
    data = data,
    outcome = "y",
    covariates = c("w1", "w2", "a")
)

# define the task for propensity score estimation
task_ps <- make_sl3_Task(
    data = data,
    outcome = "a",
    covariates = c("w1", "w2")
)

# define the task for squared outcome regression
task_outcomesq <- make_sl3_Task(
    data = data,
    outcome = "y_sq",
    covariates = c("w1", "w2", "a")
)


# train the superlearners
fit_ps <- sl_ps$train(task_ps)
fit_reg <- sl_reg$train(task_reg)
fit_regsq <- sl_reg$train(task_outcomesq)

# check the CV loss 
cv_loss_ps <- fit_ps$cv_risk(eval_fun = loss_loglik_binomial)
cv_loss_reg <- fit_reg$cv_risk(eval_fun = loss_squared_error)
cv_loss_regsq <- fit_regsq$cv_risk(eval_fun = loss_squared_error)

# check the performance

