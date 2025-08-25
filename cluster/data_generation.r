# two parameters are needed to generate the data
# 1. sample size
# 2. variability type (1 - small, 2 - medium, 3 - large)


# define the data generating function
generate_data <- function(n, effect_size_type) {
    w1 <- rnorm(n, 3, 1)
    w2 <- rbinom(n, 1, 0.5)
    propensity_score <- 1 / (exp(-(1 + 2.7 * w2 - 0.5 * w1 - 0.2 * w1 * w2)) + 1)
    a <- rbinom(n, 1, propensity_score)
    if (effect_size_type == 2) {
        sd_y <- 0.5 + 0.5 * a +
          (0.3 + 0.3 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 1) {
        sd_y <- 0.5 + a + (0.03 + 2 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 3) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.4 * a) * w1^2 + (0.12 + 0.4 * a) * w2
  } else if (effect_size_type == 4) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.03 * a) * w1^2 + (0.12 + 0.08 * a) * w2
  }
    y <- rnorm(n, 1.3 * w2 - 0.5 * log(w1^2) + 2 * a + 3.5 * a * w2, sd_y)
    y_sq <- y^2
    return(data.frame(w1, w2, a, y, y_sq))
}

# rct version
generate_data_rct <- function(n, effect_size_type) {
    w1 <- rnorm(n, 3, 1)
    w2 <- rbinom(n, 1, 0.5)
    a <- rbinom(n, 1, 0.5)
    if (effect_size_type == 2) {
        sd_y <- 0.5 + 0.5 * a +
          (0.3 + 0.3 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 1) {
        sd_y <- 0.5 + a + (0.03 + 2 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 3) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.4 * a) * w1^2 + (0.12 + 0.4 * a) * w2
  } else if (effect_size_type == 4) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.03 * a) * w1^2 + (0.12 + 0.08 * a) * w2
  }
    y <- rnorm(n, 1.3 * w2 - 0.5 * log(w1^2) + 2 * a + 3.5 * a * w2, sd_y)
    y_sq <- y^2
    return(data.frame(w1, w2, a, y, y_sq))
}


# rct with different pi_0 and pi_1
generate_data_rct_pi <- function(n, effect_size_type) {
    w1 <- rnorm(n, 3, 1)
    w2 <- rbinom(n, 1, 0.5)
    a <- rbinom(n, 1, 0.3)
    if (effect_size_type == 2) {
        sd_y <- 0.5 + 0.5 * a +
          (0.3 + 0.3 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 1) {
        sd_y <- 0.5 + a + (0.03 + 2 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 3) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.4 * a) * w1^2 + (0.12 + 0.4 * a) * w2
  } else if (effect_size_type == 4) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.03 * a) * w1^2 + (0.12 + 0.08 * a) * w2
  }
    y <- rnorm(n, 1.3 * w2 - 0.5 * log(w1^2) + 2 * a + 3.5 * a * w2, sd_y)
    y_sq <- y^2
    return(data.frame(w1, w2, a, y, y_sq))
}

# define functions to calculate the true causal effect size

# corresponding to effect size type 2
calc_causal_effect_size_medium_effect <- function(data) {
    # g-formula for the expected potential outcomes
    e_y_1 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)
    e_y_0 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)
    e_y_1_sq_cond <- (0.5 + 0.5 * 1 +
          (0.3 + 0.3 * 1) * data$w1^2 + (1.2 - 1) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)^2
    e_y_0_sq_cond <- (0.5 + 0.5 * 0 +
          (0.3 + 0.3 * 0) * data$w1^2 + (1.2 - 0) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)^2
    v_y_1 <- mean(e_y_1_sq_cond) - e_y_1^2
    v_y_0 <- mean(e_y_0_sq_cond) - e_y_0^2
    g_1 <- sum(data$a) / nrow(data)
    g_0 <- 1 - g_1
    es <- (e_y_1 - e_y_0) / sqrt(g_0 * v_y_0 + g_1 * v_y_1)
    return(es)
}

# corresponding to effect size type 1
calc_causal_effect_size_small_effect <- function(data) {
    # g-formula for the expected potential outcomes
    e_y_1 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)
    e_y_0 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)
    e_y_1_sq_cond <- (0.5 + 1 +
          (0.03 + 2 * 1) * data$w1^2 + (1.2 - 1) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)^2
    e_y_0_sq_cond <- (0.5 + 0 +
          (0.03 + 2 * 0) * data$w1^2 + (1.2 - 0) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)^2
    v_y_1 <- mean(e_y_1_sq_cond) - e_y_1^2
    v_y_0 <- mean(e_y_0_sq_cond) - e_y_0^2
    g_1 <- sum(data$a) / nrow(data)
    g_0 <- 1 - g_1
    es <- (e_y_1 - e_y_0) / sqrt(g_0 * v_y_0 + g_1 * v_y_1)
    return(es)
}

# corresponding to effect size type 3
calc_causal_effect_size_large_effect <- function(data) {
    # g-formula for the expected potential outcomes
    e_y_1 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)
    e_y_0 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)
    e_y_1_sq_cond <- (0.39 - 0.21 * 1 +
          (0.03 + 0.4 * 1) * data$w1^2 + (0.12 + 0.4 * 1) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)^2
    e_y_0_sq_cond <- (0.39 - 0.21 * 0 +
          (0.03 + 0.4 * 0) * data$w1^2 + (0.12 + 0.4 * 0) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)^2
    v_y_1 <- mean(e_y_1_sq_cond) - e_y_1^2
    v_y_0 <- mean(e_y_0_sq_cond) - e_y_0^2
    g_1 <- sum(data$a) / nrow(data)
    g_0 <- 1 - g_1
    es <- (e_y_1 - e_y_0) / sqrt(g_0 * v_y_0 + g_1 * v_y_1)
    return(es)
}

# corresponding to effect size type 4
calc_causal_effect_size_super_large_effect <- function(data) {
    # g-formula for the expected potential outcomes
    e_y_1 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)
    e_y_0 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)
    e_y_1_sq_cond <- (0.39 - 0.21 * 1 +
          (0.03 + 0.03 * 1) * data$w1^2 + (0.12 + 0.08 * 1) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)^2
    e_y_0_sq_cond <- (0.39 - 0.21 * 0 +
          (0.03 + 0.03 * 0) * data$w1^2 + (0.12 + 0.08 * 0) * data$w2)^2 +
          (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)^2
    v_y_1 <- mean(e_y_1_sq_cond) - e_y_1^2
    v_y_0 <- mean(e_y_0_sq_cond) - e_y_0^2
    g_1 <- sum(data$a) / nrow(data)
    g_0 <- 1 - g_1
    es <- (e_y_1 - e_y_0) / sqrt(g_0 * v_y_0 + g_1 * v_y_1)
    return(es)
}

# # define a function to calculate the true nuisance parameters
# calc_nuisance_params_smallv <- function(data) {
# # g-formula for the expected potential outcomes
# e_y_1 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)
# e_y_0 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)
# e_y_1_sq_cond <- (0.5 + 0.5 * 1 +
#       (0.3 + 0.3 * 1) * data$w1^2 + (1.2 - 1) * data$w2)^2 +
#       (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)^2
# e_y_0_sq_cond <- (0.5 + 0.5 * 0 +
#       (0.3 + 0.3 * 0) * data$w1^2 + (1.2 - 0) * data$w2)^2 +
#       (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)^2
# v_y_1 <- mean(e_y_1_sq_cond) - e_y_1^2
# v_y_0 <- mean(e_y_0_sq_cond) - e_y_0^2
# g_1 <- sum(data$a) / nrow(data)
# g_0 <- 1 - g_1
# es <- (e_y_1 - e_y_0) / sqrt(v_y_0) # nolint: indentation_linter.
# return(es)
# }

# calc_nuisance_params_mediumv <- function(data) {
# # g-formula for the expected potential outcomes
# e_y_1 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)
# e_y_0 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)
# e_y_1_sq_cond <- (0.5 - 0.3 * 1 +
#       (0.03 + 0.03 * 1) * data$w1^2 + (1.2 - 1) * data$w2)^2 +
#       (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)^2
# e_y_0_sq_cond <- (0.5 - 0.3 * 0 +
#       (0.03 + 0.03 * 0) * data$w1^2 + (1.2 - 0) * data$w2)^2 +
#       (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)^2
# v_y_1 <- mean(e_y_1_sq_cond) - e_y_1^2
# v_y_0 <- mean(e_y_0_sq_cond) - e_y_0^2
# g_1 <- sum(data$a) / nrow(data)
# g_0 <- 1 - g_1
# es <- (e_y_1 - e_y_0) / sqrt(v_y_0) # nolint: indentation_linter.
# return(es)
# }

# calc_nuisance_params_largev <- function(data) {
# # g-formula for the expected potential outcomes
# e_y_1 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)
# e_y_0 <- mean(1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)
# e_y_1_sq_cond <- (0.39 - 0.21 * 1 +
#       (0.03 + 0.03 * 1) * data$w1^2 + (0.12 + 0.08 * 1) * data$w2)^2 +
#       (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 1 + 3.5 * 1 * data$w2)^2
# e_y_0_sq_cond <- (0.39 - 0.21 * 0 +
#       (0.03 + 0.03 * 0) * data$w1^2 + (0.12 + 0.08 * 0) * data$w2)^2 +
#       (1.3 * data$w2 - 0.5 * log(data$w1^2) + 2 * 0 + 3.5 * 0 * data$w2)^2
# v_y_1 <- mean(e_y_1_sq_cond) - e_y_1^2
# v_y_0 <- mean(e_y_0_sq_cond) - e_y_0^2
# g_1 <- sum(data$a) / nrow(data)
# g_0 <- 1 - g_1
# es <- (e_y_1 - e_y_0) / sqrt(v_y_0) # nolint: indentation_linter.
# return(es)
# }