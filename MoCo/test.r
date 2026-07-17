generate_data <- function(n) {
    x1 <- rnorm(n, 3, 1)
    x2 <- rbinom(n, 1, 0.5)
    propensity_score <- 1 / (exp(-(1 + 2.7 * x2 - 0.5 * x1 - 0.2 * x1 * x2)) + 1)
    a <- rbinom(n, 1, propensity_score)
    m <- rnorm(n, 1 + a, 1)
    m <- exp(m)
    sd_y <- 0.5 + 0.5 * a +
        (0.3 + 0.3 * a) * x1^2 + (1.2 - a) * x2
   
    y <- rnorm(n, 1.3 * x2 - 0.5 * log(x1^2) + 2 * a + 3.5 * a * x2, sd_y)
    y_sq <- y^2
    return(data.frame(x1, x2, a, y, y_sq, m))
}

data <- generate_data(1000)
data$delta_m <- as.numeric(data$m < exp(1))
data$delta_y <- 1
library(MoCo)
result = moco(
  X = data$x1,
  A = data$a,
  M = data$m,
  Z = NULL,
  Y = data$y,
  Delta_M = data$delta_m,
  Delta_Y = data$delta_y,
  SL_library = c("SL.mean", "SL.glm", "SL.glm.interaction"),
  glm_formula = list(pMX = ".",
                     pMXZ = "."),
  cross_fit = TRUE,
  HAL_pMX = FALSE,
  HAL_pMXZ = FALSE,
  cv_folds = 5,
  seed_rgn = 1,
  test = TRUE,
  fwer = 0.05
)
