# this script is used to check the performance of the ranger learner
library(ranger)
library(sl3)


source("data_generation.r")

set.seed(1234)
data1 <- generate_data(n = 100, effect_size_type = 3)
data2 <- generate_data(n = 500, effect_size_type = 3)
data3 <- generate_data(n = 1000, effect_size_type = 3)
data4 <- generate_data(n = 2000, effect_size_type = 3)
data5 <- generate_data(n = 5000, effect_size_type = 3)
data6 <- generate_data(n = 10000, effect_size_type = 3)

# define a function to run the estimation using ranger
# store bias for nuisance parameters
estimate_ranger <- function(data) {
    # define the learner
    lrnr_ranger <- Lrnr_ranger$new()

    # define the task for propensity score estimation
    task_ps <- make_sl3_Task(
        data = data,
        outcome = "a",
        covariates = c("w1", "w2")
    )
    
    # define the task for outcome regression
    task_reg <- make_sl3_Task(
        data = data,
        outcome = "y",
        covariates = c("w1", "w2", "a")
    )

    # define the task for outcome regression with squared outcome
    task_regsq <- make_sl3_Task(
        data = data,
        outcome = "y_sq",
        covariates = c("w1", "w2", "a")
    )
    
    # train the learners
    fit_ps <- lrnr_ranger$train(task_ps)
    fit_reg <- lrnr_ranger$train(task_reg)
    fit_regsq <- lrnr_ranger$train(task_regsq)

    # get the bias for the nuisance parameters
    bias_ps <- mean(fit_ps$predict(task_ps) - data$a)
    bias_reg <- mean(fit_reg$predict(task_reg) - data$y)
    bias_regsq <- mean(fit_regsq$predict(task_regsq) - data$y_sq)
    
    # return the bias
    return(list(bias_ps = bias_ps, bias_reg = bias_reg, bias_regsq = bias_regsq))
}

# run the estimation
bias_ps1 <- estimate_ranger(data1)$bias_ps
bias_reg1 <- estimate_ranger(data1)$bias_reg
bias_regsq1 <- estimate_ranger(data1)$bias_regsq

bias_ps2 <- estimate_ranger(data2)$bias_ps
bias_reg2 <- estimate_ranger(data2)$bias_reg
bias_regsq2 <- estimate_ranger(data2)$bias_regsq

bias_ps3 <- estimate_ranger(data3)$bias_ps
bias_reg3 <- estimate_ranger(data3)$bias_reg
bias_regsq3 <- estimate_ranger(data3)$bias_regsq

bias_ps4 <- estimate_ranger(data4)$bias_ps
bias_reg4 <- estimate_ranger(data4)$bias_reg
bias_regsq4 <- estimate_ranger(data4)$bias_regsq


bias_ps5 <- estimate_ranger(data5)$bias_ps
bias_reg5 <- estimate_ranger(data5)$bias_reg
bias_regsq5 <- estimate_ranger(data5)$bias_regsq

bias_ps6 <- estimate_ranger(data6)$bias_ps
bias_reg6 <- estimate_ranger(data6)$bias_reg
bias_regsq6 <- estimate_ranger(data6)$bias_regsq

# output the bias in the terminal in a table
bias_table <- data.frame(
    n = c(100, 500, 1000, 2000, 5000, 10000),
    bias_ps = c(bias_ps1, bias_ps2, bias_ps3, bias_ps4, bias_ps5, bias_ps6),
    bias_reg = c(bias_reg1, bias_reg2, bias_reg3, bias_reg4, bias_reg5, bias_reg6),
    bias_regsq = c(bias_regsq1, bias_regsq2, bias_regsq3, bias_regsq4, bias_regsq5, bias_regsq6)
)
print(bias_table)





