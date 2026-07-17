# define the data generating function
# four types of effect sizes:
# 1. small
# 2. medium
# 3. large
# 4. huge

# three types of propensity scores:
# 1. no confounding with p=0.5
# 2. no confounding with p=0.3
# 3. confounding

generate_data <- function(n, effect_size_type, propensity_score_type) {
    w1 <- rnorm(n, 3, 1)
    w2 <- rbinom(n, 1, 0.5)
    if (propensity_score_type == 1) {
        propensity_score <- 0.5
    } else if (propensity_score_type == 2) {
        propensity_score <- 0.3
    } else if (propensity_score_type == 3) {
        propensity_score <- 1 / (exp(-(1 + 2.7 * w2 - 0.5 * w1 - 0.2 * w1 * w2)) + 1)
    } else {
        stop("Invalid propensity score type")
    }
    a <- rbinom(n, 1, propensity_score)
    if (effect_size_type == 1) {
        sd_y <- 0.5 + a + (0.03 + 2 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 2) {
        sd_y <- 0.5 + 0.5 * a +
          (0.3 + 0.3 * a) * w1^2 + (1.2 - a) * w2
    } else if (effect_size_type == 3) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.4 * a) * w1^2 + (0.12 + 0.4 * a) * w2
  } else if (effect_size_type == 4) {
    sd_y <- 0.39 - 0.21 * a +
      (0.03 + 0.03 * a) * w1^2 + (0.12 + 0.08 * a) * w2
  } else {
    stop("Invalid effect size type")
  }
    y <- rnorm(n, 1.3 * w2 - 0.5 * log(w1^2) + 2 * a + 3.5 * a * w2, sd_y)
    y_sq <- y^2
    return(data.frame(w1, w2, a, y, y_sq))
}